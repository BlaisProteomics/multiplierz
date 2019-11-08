# Copyright 2008 Dana-Farber Cancer Institute
# multiplierz is distributed under the terms of the GNU Lesser General Public License
#
# This file is part of multiplierz.
#
# multiplierz is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# multiplierz is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with multiplierz.  If not, see <http://www.gnu.org/licenses/>.

import os
import re
import wx

#import wxmpl

from numpy import array, hypot

MZ_EXT = ('.raw', '.wiff', '.mzml', '.mzml.gz')
MZ_EXT_2 = MZ_EXT + tuple((e + '.lnk') for e in MZ_EXT) # with shortcuts included

MZ_WILDCARD = 'MS Data Files (%s)|%s' % ('; '.join(('*' + e) for e in MZ_EXT),
                                         '; '.join(('*' + e) for e in MZ_EXT))


class mzApp():
    def __init__(self):
        #self.app = wx.PySimpleApp()
        self.app = wx.App(False)

    def launch(self):
        self.app.MainLoop()


#class mzForm(wx.Frame):
    #def __init__(self, parent=None, title="mzForm", items=None, function=None, size=(600,450)):
        #wx.Frame.__init__(self, parent, -1, title, size=size)

        #self.SetIcon(wx.Icon(os.path.join(install_dir, 'images', 'icons', 'multiplierz.ico'),
                             #wx.BITMAP_TYPE_ICO))

        #items = items or []

        #self.function = function
        #self.files = dict()
        #self.variables = set()

        #self.pane = wx.Panel(self, -1)

        ## status bar
        #self.statusbar = self.CreateStatusBar()
        #self.statusbar.SetFieldsCount(2)
        #self.statusbar.SetStatusText("Ready", 0)

        #gbs = wx.GridBagSizer(10,15)

        #for i,(ctrl_type, label, variable, value) in enumerate(items):
            #self.variables.add(variable)

            #gbs.Add( wx.StaticText(self.pane, -1, label, style=wx.ALIGN_RIGHT),
                     #(i,0), flag=wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT )

            #ctrl_type = ctrl_type.lower()

            #if ctrl_type == 'file' or ctrl_type == 'files':
                #txt_ctrl = wx.TextCtrl(self.pane, -1, value, name=variable)
                #gbs.Add( txt_ctrl,
                         #(i,1), flag=wx.EXPAND )

                #browse_btn = wx.Button(self.pane, -1, 'Browse')
                #gbs.Add( browse_btn,
                         #(i,2), flag=wx.EXPAND )

                #if ctrl_type == 'file':
                    #browse_btn.Bind(wx.EVT_BUTTON, lambda e: self.on_file_browse(txt_ctrl))
                    #self.files[str(variable)] = ''
                #else:
                    #browse_btn.Bind(wx.EVT_BUTTON, lambda e: self.on_files_browse(txt_ctrl))
                    #self.files[str(variable)] = []
            #else:
                #if ctrl_type == 'text':
                    #gbs.Add( wx.TextCtrl(self.pane, -1, value, name=variable),
                             #(i,1), (1,2), flag=wx.EXPAND|wx.ALIGN_CENTER )
                #elif ctrl_type == 'radio':
                    #gbs.Add( wx.RadioBox(self.pane, -1, label='', choices=value, name=variable),
                             #(i,1), (1,2), flag=wx.EXPAND|wx.ALIGN_CENTER )
                #elif ctrl_type == 'check':
                    #check_box = wx.CheckBox(self.pane, -1, label='', name=variable)
                    #gbs.Add( check_box,
                             #(i,1), (1,2), flag=wx.EXPAND|wx.ALIGN_CENTER )
                    #check_box.SetValue(value)

        ## submit button
        #submit_btn = wx.Button(self.pane, -1, "Submit")
        #gbs.Add( submit_btn,
                 #(i+1,0), (1,3), flag=wx.ALIGN_CENTER )
        #self.Bind(wx.EVT_BUTTON, self.on_submit, submit_btn)

        #gbs.AddGrowableCol(1, 1)
        #gbs.AddGrowableRow(i+1, 1)

        #box = wx.BoxSizer()
        #box.Add(gbs, 1, wx.ALL|wx.EXPAND, 5)
        #self.pane.SetSizerAndFit(box)

    #def on_file_browse(self, text_ctrl):
        #file_chooser = wx.FileDialog(None, "Choose File:", style=wx.FD_OPEN)

        #if file_chooser.ShowModal() == wx.ID_OK:
            #self.files[str(text_ctrl.GetName())] = file_chooser.GetPath()
            #text_ctrl.SetValue(file_chooser.GetPath())
            #text_ctrl.SetToolTipString(file_chooser.GetPath())
        #else:
            #self.files[str(text_ctrl.GetName())] = ''
            #text_ctrl.SetValue('')
            #text_ctrl.SetToolTipString('')

        #file_chooser.Destroy()

    #def on_files_browse(self, text_ctrl):
        #file_chooser = wx.FileDialog(None, "Choose Files:", style=wx.FD_MULTIPLE)

        #if file_chooser.ShowModal() == wx.ID_OK:
            #self.files[str(text_ctrl.GetName())] = file_chooser.GetPaths()
            #text_ctrl.SetValue('; '.join(file_chooser.GetPaths()))
            #text_ctrl.SetToolTipString('\n'.join(file_chooser.GetPaths()))
        #else:
            #self.files[str(text_ctrl.GetName())] = []
            #text_ctrl.SetValue('')
            #text_ctrl.SetToolTipString('')

        #file_chooser.Destroy()

    #def on_submit(self, event):
        #self.statusbar.SetStatusText("Working...", 0)
        #wx.BeginBusyCursor(wx.HOURGLASS_CURSOR)

        ##Get Variables:
        #values = self.files.copy()

        #for c in self.pane.GetChildren():
            #if c.GetName() in self.variables and c.GetName() not in self.files:
                #if isinstance(c, wx.RadioBox):
                    #values[str(c.GetName())] = c.GetStringSelection()
                #elif not isinstance(c, wx.StaticText):
                    #values[str(c.GetName())] = c.GetValue()

        #self.function(**values)
        #self.statusbar.SetStatusText("Ready", 0)
        #wx.EndBusyCursor()


#class mzPlot(wx.Frame):
    #def __init__(self, parent=None, title="mzPlot", size=(600,450)):
        #wx.Frame.__init__(self, parent, -1, title, size=size)

        #self.xy_data = None
        #self.last_anno = None
        #self.tooltip_str = '%%3.1f, %%3d' # default tooltip string

        ##Icon
        #self.SetIcon(wx.Icon(os.path.join(install_dir, 'images', 'icons', 'multiplierz.ico'),
                             #wx.BITMAP_TYPE_ICO))

        ##add menu bar
        #menu_bar = wx.MenuBar()

        ##Edit Menu
        #edit_menu = wx.Menu()

        #change_title = edit_menu.Append(-1, 'Change &Title\tCtrl+T', 'Change Plot Title')
        #self.Bind(wx.EVT_MENU, self.on_title, change_title)

        #x_label = edit_menu.Append(-1, 'Change &X Axis Label\tCtrl+X', 'Change X Axis Label')
        #self.Bind(wx.EVT_MENU, self.on_xlabel, x_label)

        #y_label = edit_menu.Append(-1, 'Change &Y Axis Label\tCtrl+Y', 'Change Y Axis Label')
        #self.Bind(wx.EVT_MENU, self.on_ylabel, y_label)

        #menu_bar.Append(edit_menu, "&Edit")

        #save_menu = wx.Menu()

        #save_image = save_menu.Append(-1, '&Save Image\tCtrl+S', 'Save Plot as Image')
        #self.Bind(wx.EVT_MENU, self.on_save, save_image)

        #menu_bar.Append(save_menu, "&Save")

        #resize_menu = wx.Menu()

        #resize_800 = resize_menu.Append(-1, "800x600\tAlt+1", "Resize Plot to 800x600")
        #self.Bind(wx.EVT_MENU, lambda e: self.on_resize((800,600)), resize_800)

        #resize_1200 = resize_menu.Append(-1, "1200x900\tAlt+2", "Resize Plot to 1200x900")
        #self.Bind(wx.EVT_MENU, lambda e: self.on_resize((1200,900)), resize_1200)

        #resize_1400 = resize_menu.Append(-1, "1400x1050\tAlt+3", "Resize Plot to 1400x1050")
        #self.Bind(wx.EVT_MENU, lambda e: self.on_resize((1400,1050)), resize_1400)

        #menu_bar.Append(resize_menu, "&Resize")

        #self.SetMenuBar(menu_bar)

        #self.plot_panel = wxmpl.PlotPanel(self, -1, (1.6, 1.2))
        #self.plot_panel.mpl_connect('button_release_event', self.on_click)

        #self.figure = self.plot_panel.get_figure()
        #a = self.figure.add_axes([0.125, 0.1, 0.775, 0.8])
        #a.set_title(title)

        #self.plot_panel.draw()

        #box = wx.BoxSizer()
        #box.Add(self.plot_panel, 1, wx.EXPAND, 0)
        #self.SetSizerAndFit(box)
        #self.SetSize(size)

    #def on_resize(self, size):
        #self.SetSize(size)
        #self.SendSizeEvent()

    #def on_title(self, event):
        #with wx.TextEntryDialog(self, 'Title this graph',
                                #'Enter Graph Title',
                                #self.GetTitle()) as title_dlg:
            #if title_dlg.ShowModal() == wx.ID_OK:
                #title = title_dlg.GetValue()
                #self.SetTitle(title)
                #self.figure.get_axes()[0].set_title(title)
                #self.plot_panel.draw()

    #def on_xlabel(self, event):
        #with wx.TextEntryDialog(self, 'Change X-Axis Label',
                                #'Enter X-Axis Label',
                                #self.figure.get_axes()[0].get_xlabel()) as xlabel_dlg:
            #if xlabel_dlg.ShowModal() == wx.ID_OK:
                #title = xlabel_dlg.GetValue()
                #self.figure.get_axes()[0].set_xlabel(title)
                #self.plot_panel.draw()

    #def on_ylabel(self, event):
        #with wx.TextEntryDialog(self, 'Change Y-Axis Label',
                                #'Enter Y-Axis Label',
                                #self.figure.get_axes()[0].get_ylabel()) as ylabel_dlg:
            #if ylabel_dlg.ShowModal() == wx.ID_OK:
                #title = ylabel_dlg.GetValue()
                #self.figure.get_axes()[0].set_ylabel(title)
                #self.plot_panel.draw()

    #def on_save(self, event):
        #wildcard = ("PNG (*.png)|*.png|"
                    #"PDF (*.pdf)|*.pdf|"
                    #"PS (*.ps)|*.ps|"
                    #"EPS (*.eps)|*.eps|"
                    #"SVG (*.svg)|*.svg")
        #formats = ('PNG', 'PDF', 'PS', 'EPS', 'SVG')

        #with wx.FileDialog(self, "Save figure as...",
                           #wildcard=wildcard, style=wx.FD_SAVE) as dlg:
            #if dlg.ShowModal() == wx.ID_OK:
                #self.plot_panel.print_figure(dlg.GetPath(),
                                             #format=formats[dlg.GetFilterIndex()])

    #def closest_point(self, event):
        #if self.xy_data is None:
            #return None

        #axes = event.canvas.figure.get_axes()[0]

        #xlim = axes.get_xlim()
        #ylim = axes.get_ylim()

        #xy_data = [(x,y) for x,y in self.xy_data
                   #if xlim[0] <= x <= xlim[1] and ylim[0] <= y <= ylim[1]]

        #if not xy_data:
            #return None

        #e_xy = array([event.x, event.y])

        #xy = min((axes.transData.transform([x,y]) for x,y in xy_data),
                 #key = lambda xy: hypot(*(e_xy - xy)))

        ## 10 pixel threshold for labeling
        #if all(abs(xy - e_xy) < 10.0):
            #return (tuple(abs(axes.transData.inverted().transform(xy))),
                    #tuple(axes.transData.inverted().transform(xy+5)))
        #else:
            #return None

    #def on_click(self, event):
        #'''Annotate the point closest to the cursor if it is within range'''

        #if event.inaxes:
            #xy_o = self.closest_point(event)
            #if xy_o:
                #xy,o = xy_o

                #if self.last_anno is not None:
                    #self.last_anno.remove()

                #tip = self.tooltip_str % xy

                #axes = self.figure.get_axes()[0]

                #t = axes.text(o[0], o[1], tip)
                #self.last_anno = t
                #event.canvas.draw()

                #return

        #if self.last_anno is not None:
            #self.last_anno.remove()
            #self.last_anno = None

        #event.canvas.draw()

    #def plot(self, *args, **kwargs):
        #'''A simple wrapper for matplotlib's axes.plot() function. If you
        #want to do something more complicated, you can access the figure
        #directly using mzPlot.figure'''

        #self.figure.clear()
        #axes = self.figure.add_axes([0.125, 0.1, 0.775, 0.8])
        #self.xy_data = axes.plot(*args, **kwargs)[0].get_xydata()

        #self.plot_panel.draw()

    #def plot_xic(self, title="XIC", data=None, scan_dot=None, other_MS2s=None):
        #if data is None:
            #raise TypeError("Required argument 'data' cannot be None")

        #self.tooltip_str = '(%%3.%df, %%3.%df)' % (settings.xic_time_figs,
                                                   #settings.xic_int_figs)

        #mz_image._make_xic(self.figure, None,
                           #[x for x,y in data],
                           #[y for x,y in data],
                           #scan_dot,
                           #[x for x,y in other_MS2s] if other_MS2s else [],
                           #[y for x,y in other_MS2s] if other_MS2s else [],
                           #title)

        #self.plot_panel.draw()

    #def plot_full_ms(self, title="Full MS", scan=None, scan_mz=None):
        #if scan is None:
            #raise TypeError("Required argument 'scan' cannot be None")

        #self.tooltip_str = '(%%3.%df, %%3.%df)' % (settings.ms1_mz_figs,
                                                   #settings.ms1_int_figs)

        #mz_image._make_ms1(self.figure,
                           #None,
                           #scan,
                           #scan.mode,
                           #[scan_mz] if scan_mz else None,
                           #title,
                           #settings.MS1_view_mz_window / 2)

        #self.plot_panel.draw()

    #def plot_ms_ms(self, title="MS-MS", scan=None):
        #if scan is None:
            #raise TypeError("Required argument 'scan' cannot be None")

        #self.tooltip_str = '(%%3.%df, %%3.%df)' % (settings.ms2_mz_figs,
                                                   #settings.ms2_int_figs)

        #mz_image._make_ms2(self.figure,
                           #scan,
                           #scan.mode,
                           #None,
                           #title=title)

        #self.plot_panel.draw()

    #def plot_venn(self, A, B, AB, A_label, B_label, title='Venn Diagram', eps=0.001):
        #'''Plot a proportional 2-set Venn diagram. A and B are the sizes of the two sets,
        #AB is the size of the intersection, and eps is an error margin for the proportional
        #placement. E.g. if eps is 0.01 then the areas of the plot will be accurate to ~1%.

        #A lower eps will give a more accurate plot at the expense of longer running time.
        #The method uses a bisecting search algorithm to find the right proportions.'''

        #mz_image.make_venn(self.figure, A, B, AB, A_label, B_label, title, eps)

        #self.plot_panel.draw()


class NumValidator(wx.PyValidator):
    '''This is a generic float-validating class which is used in various
    forms to require a non-negative floating-point number from a text control.
    A flag indicates whether 0.0 is considered valid.

    If the validator fails, it will highlight the control and pop up an
    error message. The control should have a reasonable name to stick
    in the message text.
    '''
    def __init__(self, func=float, flag=False):
        wx.PyValidator.__init__(self)
        self.flag = not flag # reverse the flag for convenience later
        self.func = func # used to cast the result, to allow int or float validation

    def Clone(self):
        return NumValidator(self.func, not self.flag)

    def Validate(self, win):
        tc = self.GetWindow()
        val = tc.GetValue()
        nm = tc.GetName()

        try:
            v = self.func(val)
            if self.flag and v <= 0.0:
                wx.MessageBox("%s must be positive" % nm, "Error")
            elif v < 0.0:
                wx.MessageBox("%s must be non-negative" % nm, "Error")
            else:
                tc.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
                tc.Refresh()
                return True
        except ValueError:
            wx.MessageBox("%s is not a valid number" % nm, "Error")

        tc.SetBackgroundColour("YELLOW")
        tc.Clear()
        tc.SetFocus()
        tc.Refresh()
        return False


class ProgressBar:
    def __init__(self, title, entries):
        self.progressBar = wx.ProgressDialog(title, "Time remaining", entries, style=wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME | wx.PD_AUTO_HIDE | wx.PD_CAN_ABORT | wx.PD_APP_MODAL)

    def update(self, count):
        self.progressBar.Update(count)

        if self.progressBar.Update(count)[0] == 0:
            cancelMsg = wx.MessageDialog(None, "Are you sure you want to cancel?",'Continue?',wx.YES_NO | wx.ICON_QUESTION)
            cancelMsgAnswer = cancelMsg.ShowModal()
            if cancelMsgAnswer == wx.ID_YES:
                return False
            else:
                self.progressBar.Resume()

        return True

    def destroy(self):
        self.progressBar.Destroy()


def alerts(message='multiplierz', title='multiplierz', method='popup'):
    """Alert system for displaying popup or writing message to a file

    Available methods equal 'popup' or 'file'.
    If a file type method is chosen, the title represents the file name and location.

    Example:
    >> alerts('popup', 'Test', 'Hello. The test worked. Please click OK.')

    """

    message = str(message)
    title = str(title)
    if method == 'popup':
        try:
            wx.MessageBox(message, title)
        except wx._core.PyNoAppError:
            app = mzApp()
            app.launch()
            wx.MessageBox(message, title)
    elif method == 'file':
        fh = open(title,'w')
        fh.write(message)
        fh.close()


def file_chooser(title='Choose a file:', default_path = '', default_file = '',
                 mode='r', wildcard='*', parent_obj = None):
    """Provides a file chooser dialog and returns the file path(s) when the file(s) is selected.
    mode option provides file dialog type: read single, read multiple, or save single.
    mode = 'r' creates an 'Open' file dialog for single file read.
    mode = 'm' creates an 'Open' file dialog for multiple files read.
    mode = 'w' creates a 'Save' file dialog.

    wildcard option can be specified as "*.xls"

    Example:
    >> file_chooser(title='Choose a file:', mode='r')

    """

    #wildcard = "%s|%s" % (wildcard, wildcard)

    style = { 'r': wx.FD_OPEN,
              'm': wx.FD_MULTIPLE,
              'w': wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT }[mode]

    try:
        file_chooser = wx.FileDialog(parent_obj, title, wildcard=wildcard, style=style,
                                     defaultDir = default_path, defaultFile = default_file)
    except wx._core.PyNoAppError:
        app = mzApp()
        app.launch()
        file_chooser = wx.FileDialog(parent_obj, title, wildcard=wildcard, style=style,
                                     defaultDir = default_path, defaultFile = default_file)


    file_name = None
    if file_chooser.ShowModal() == wx.ID_OK:
        if mode == 'm':
            file_name = file_chooser.GetPaths()
        else:
            file_name = file_chooser.GetPath()
    file_chooser.Destroy()

    return file_name

def report_chooser(title=None, mode='r', parent = None, **kwargs):
    '''A specialized file_chooser function for multiplierz files. Otherwise,
    works just like file_chooser.

    'parent' is the parent of the dialog--used when this is called by a GUI
    element. It is perfectly fine to leave it as None, and the GUI frame will
    have no parent.

    'title' can be left blank for the following default titles based on mode:
    r - 'Choose multiplierz File:'
    w - 'Save File:'
    m - 'Choose multiplierz Files:'

    'mode' is one of 'r', 'w', and 'm', just as for file_chooser.

    **kwargs can include any additional options to pass to the FileDialog constructor,
    such as defaultDir (default directory).'''

    # For legacy reasons, these are often misplaced in scripts.
    # But they're both necessarily typed differently, so its sortable.
    if isinstance(parent, str) and not isinstance(title, str):
        title, parent = parent, title

    wildcard = ("Worksheets (*.xls; *.xlsx)|*.xls; *.xlsx|"
                "Comma-separated Values (*.csv)|*.csv|"
                "mzResults Database (*.mzd)|*.mzd|"
                "mzIdentML (*.mzid)|*.mzid")

    if not title:
        title = {'r': 'Choose multiplierz File:',
                 'w': 'Save File:',
                 'm': 'Choose multiplierz Files:'}[mode]

    style = { 'r': wx.FD_OPEN,
              'm': wx.FD_MULTIPLE,
              'w': wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT }[mode]

    #index = {'.xls': 0,
             #'.xlsx': 0,
             #'.csv': 1,
             #'.mzd': 2}[settings.default_format]

    index = 0

    try:
        file_dialog = wx.FileDialog(parent, title, wildcard=wildcard, style=style, **kwargs)
    except wx._core.PyNoAppError as err:
        app = mzApp()
        app.launch()
        file_dialog = wx.FileDialog(None, title, wildcard=wildcard, style=style, **kwargs)
    file_dialog.SetFilterIndex(index)

    file_name = None
    if file_dialog.ShowModal() == wx.ID_OK:
        if mode == 'm':
            file_name = file_dialog.GetPaths()
        else:
            file_name = file_dialog.GetPath()
    file_dialog.Destroy()

    return file_name




class FileArrayDialog(wx.Dialog):
    def __init__(self, parent, filetypes, fileLimit = None):
        wx.Dialog.__init__(self, parent, -1, title = 'Select Input Files',
                           style = wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        
        self.filetypes = {}
        
        self.fileArray = wx.ListCtrl(self, -1, style = wx.LC_REPORT | wx.LC_EDIT_LABELS | wx.LC_HRULES | wx.LC_VRULES)
        for i, (filetype, ext) in enumerate(filetypes):
            self.fileArray.InsertColumn(i, filetype)

            assert filetype not in self.filetypes, "Non-unique file identifiers! %s" % filetype
            self.filetypes[filetype] = ext
            
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.getFile, self.fileArray)
    
        self.goButton = wx.Button(self, -1, "Use Selected Files")
        
        self.Bind(wx.EVT_BUTTON, self.complete, self.goButton)
        self.Bind(wx.EVT_CLOSE, self.abort)
        
        
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.fileArray, 1, wx.ALL | wx.EXPAND, 10)
        box.Add(self.goButton, 0, wx.ALL | wx.EXPAND, 10)
        self.SetSizerAndFit(box)
        
        
        self.SetSize(wx.Size(1200, 300))
        self.Bind(wx.EVT_SIZE, self.resize)

        for i in range(0, 10):
            self.fileArray.Append([''] * self.fileArray.GetColumnCount())
        self.resize(None)
        self.CentreOnScreen()
    
    def resize(self, event):
        arraywidth = self.fileArray.GetSize()[0] - 10
        for col in range(self.fileArray.GetColumnCount()):
            self.fileArray.SetColumnWidth(col, arraywidth / self.fileArray.GetColumnCount()) 
        if event:
            event.Skip()
            
    def getFile(self, event):
        row = event.GetIndex()
        for col in range(self.fileArray.GetColumnCount()):
            filetype = self.fileArray.GetColumn(col).GetText()
            exts = self.filetypes[filetype]
            wildcard = '|'.join(['%s|%s' % (x, x) for x in exts] + ['*|*']) 
            #wildcard = '|'.join(exts + ['*'])
            givenfile = file_chooser(title = 'Choose %s:' % filetype, wildcard = wildcard)
            self.fileArray.SetStringItem(index = row, col = col, label = givenfile)
            
        
    def returnFiles(self):
        filesets = []
        for r in range(self.fileArray.GetItemCount()):
            fileset = []
            for c in range(self.fileArray.GetColumnCount()):
                fileset.append(self.fileArray.GetItem(r, c).GetText())
            if any(fileset):
                filesets.append(fileset)
        
        return filesets
    
    def complete(self, event):
        self.EndModal(wx.ID_OK)
    def abort(self, event):
        raise RuntimeError("User cancelled file selection.")
    
def open_filearray(parent = None, filetypes = None):
    assert filetypes and all([len(x) == 2 for x in filetypes])
    if parent:
        dialog = FileArrayDialog(parent, filetypes)
    else:
        app = mzApp()
        app.launch()
        dialog = FileArrayDialog(None, filetypes)
    dialog.ShowModal()
    files = dialog.returnFiles()
    return files
        
        

def text_input(prompt='', value='', title=''):
    """Provides a text input dialog and returns the user input

    The value field can be used to enter the default initial value.

    Example:
    >> textInput('Enter your name:', title='Name')

    """

    try:
        dlg = wx.TextEntryDialog(None, message=prompt, caption=title, defaultValue=value)
    except wx._core.PyNoAppError:
        app = mzApp()
        app.launch()
        dlg = wx.TextEntryDialog(None, message=prompt, caption=title, defaultValue=value)

    dlg.ShowModal()
    output = dlg.GetValue()
    dlg.Destroy()

    return output






def directory_chooser(parent = None, title = None):
    try:
        dialog = wx.DirDialog(None)
    except wx._core.PyNoAppError:
        app = mzApp()
        app.launch()
        dialog = wx.DirDialog(None)
        
    if dialog.ShowModal() == wx.ID_OK:
        return dialog.GetPath()
    else:
        return None
