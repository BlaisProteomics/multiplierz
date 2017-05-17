from multiplierz.mzAPI import mzFile
import sys



# Its not clear if sys.getsizeof() will work on mzScan objects or etc; careful!



# Substantial voodoo; the class stores cached data in the scope of a
# multiprocess-threaded function; the function loops, responding to commands
# or otherwise loading data unless it has hit the size cap (or exhausted
# the data file of easily cachable things.)

# For a proof-of-concept, at least, don't do the autonomous-data-intake;
# getting results from calls and caching results from preorder commands is
# enough.
#
# It may be tricky making those interleave properly; a call can wait while
# preorders or whatever else retrieve data, but the next output on the queue
# has to always be from the most recent command (we *presume* no ones going to
# be silly enough to hand two copies of the main class to different threads.)


# The thread should only block on the command queue when there's nothing else
# left for it to do!


def arghash(method, args):
    """Argh!"""
    foo, bar = args
    bar = tuple((x, y) for x, y in bar.items())
    return hash((method, foo, bar))

def async_mzFile_internal(datafile, size_cap, input, output):
    try:
        from collections import deque
        from Queue import Empty
        
        data = mzFile(datafile)
        cacheRecord = deque()
        cache = {}
        commands = []
        
        while True:
            if not commands:
                commands = [input.get(block = True)]
                
            while True:
                try:
                    commands.append(input.get_nowait())
                except Empty:
                    break
            
            if any(x[0] == 'close' for x in commands):
                break
                
            
            call = next((x for x in commands if x[0] == 'call'), None)
            if call:
                commands.remove(call)

                _, method, argkwarg = call
                callhash = arghash(method, argkwarg)
                if callhash in cache:
                    output.put((call, cache[callhash]))
                else:
                    args, kwargs = argkwarg
                    returnval = getattr(data, method)(*args, **kwargs)
                    cache[callhash] = returnval
                    output.put((call, returnval))
            
            else:
                assert commands
                com = commands[0]
                _, method, argkwarg = com
                callhash = arghash(method, argkwarg)
                commands = commands[1:]
                if callhash in cache:
                    continue
                args, kwargs = argkwarg
                returnval = getattr(data, method)(*args, **kwargs)
                cache[callhash] = returnval
        
        return
    except Exception as err:
        import traceback
        traceback.print_exc()
        print '------------------'
        raise err        
            
            
        
        

class async_mzFile(object):
    def __init__(self, filename, size_cap = 1000 * 1000 * 1000):
        import multiprocessing
        
        self.data_file = filename
        
        self.commands = multiprocessing.Queue()
        self.results = multiprocessing.Queue()
        
        self.source = multiprocessing.Process(target = async_mzFile_internal,
                                              args = (filename, size_cap,
                                                      self.commands, self.results))
        self.source.start()
        
        
    def proxycall(self, method, *args, **kwargs):
        if not self.source.is_alive():
            err, returnval = self.source.get()
            print "Data retrieval thread crashed with the following error:\n\n"
            print returnval
            raise err
            
        
        allargs = args, kwargs
        
        self.commands.put(('call', method, allargs))
        
        foo = self.results.get()
        command, returnval = foo
        assert command == ('call', method, allargs)
        
        return returnval
    
    def preorder(self, method, *args, **kwargs):
        if not self.source.is_alive():
            err, returnval = self.source.get()
            print "Data retrieval thread crashed with the following error:\n\n"
            print returnval
            raise err        
        
        allargs = args, kwargs
        
        self.commands.put(('preorder', method, allargs))
        #ack = self.results.get()
        #assert ack == (method, allargs)
        
    def __getattr__(self, method):
        # Non-method attributes may be problematic.
        
        #if not self.check_method(method):
            #raise NotImplementedError, "Only method attributes are represented in async_mzFile objects!"
        
        def proxymethod(*args, **kwargs):
            return self.proxycall(method, *args, **kwargs)
        return proxymethod
    

    def close(self):
        self.commands.put(('close', None, None))
        #self.source.wait()
        self.source.join()
        
    def __del__(self):
        if self.source.is_alive():
            self.close()
        
        
        
        
        
        
        
        
if __name__ == '__main__':
    from multiplierz.mzAPI.mzMemo import async_mzFile
    from time import clock, sleep
    foo = async_mzFile('example_file.raw')
    start = clock()
    foo.preorder('xic', 0, 999999)
    print clock() - start
    
    sleep(3)
    start = clock()
    bar = foo.xic(0, 999999)
    print clock() - start
    
    foo.close()