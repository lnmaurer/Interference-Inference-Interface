import Tkinter, tkFileDialog
import ttk
import tkMessageBox

class Interface:
  """The class for the GUI interface"""
  
  viewWidth = 500 #width of the view canvas
  viewHeight = 100 #height of the view canvas

  def __init__(self):  
    self.root = Tkinter.Tk()
    self.root.title("Leon's Olde Interference & Diffraction Simulator")
    
    #The menubar and menus
    menubar = Tkinter.Menu(self.root)
  
    #the file menu
    filemenu = Tkinter.Menu(menubar, tearoff=0)
    filemenu.add_command(label="Save Experiment As", accelerator="Ctrl+S")
    filemenu.add_command(label="Load Experiment", accelerator="Ctrl+O")
    filemenu.add_separator()
    filemenu.add_command(label="Exit", accelerator="Ctrl+Q", command=self.root.quit)
    #bind keys to the actions
    #self.root.bind_all('<Control-s>', lambda arg: self.saveExperiment())
    #self.root.bind_all('<Control-o>', lambda arg: self.loadExperiment())
    self.root.bind_all('<Control-q>', lambda arg: self.root.quit())
    menubar.add_cascade(label="File", menu=filemenu)
  
    #the edit menu
    editmenu = Tkinter.Menu(menubar, tearoff=0)
    editmenu.add_command(label="Cut", accelerator="Ctrl+X") #todo: add command???
    editmenu.add_command(label="Copy", accelerator="Ctrl+C") #todo: add command???
    editmenu.add_command(label="Paste", accelerator="Ctrl+V") #todo: add command???
    menubar.add_cascade(label="Edit", menu=editmenu)
  
    self.root.config(menu=menubar)    
    
    
if __name__ == "__main__":
#the root window
  gui = Interface() #make the interface
  gui.root.mainloop() #set everything in motion