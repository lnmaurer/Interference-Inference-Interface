import Tkinter, tkFileDialog
import ttk
import tkMessageBox
import numpy as np
from PIL import Image, ImageTk

class Interface:
  """The class for the GUI interface"""
  
  viewWidth = 800 #width of the view canvas
  viewHeight = 400 #height of the view canvas
  minFreq = 1
  maxFreq = 37
  barrierX = 100

  def __init__(self):    
#The root window
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
    editmenu.add_command(label="Cut", accelerator="Ctrl+X") #todo: add command
    editmenu.add_command(label="Copy", accelerator="Ctrl+C") #todo: add command
    editmenu.add_command(label="Paste", accelerator="Ctrl+V") #todo: add command
    menubar.add_cascade(label="Edit", menu=editmenu)
  
    self.root.config(menu=menubar)
    
#The view frame
    self.viewFrame = ttk.Labelframe(self.root, text='View')
    self.viewFrame.grid(column=0,row=0,sticky='nsew',padx=5,pady=5)
    
    self.canvas = Tkinter.Canvas(self.viewFrame, width=self.viewWidth, height=self.viewHeight)
    self.canvas.grid(column=1, row=0, columnspan=3, rowspan=3, sticky='nsew', padx=5, pady=5)    
    
#The barrier frame frame and intial barrier setup
    self.barrierFrame = ttk.Labelframe(self.root, text='Barrier')
    self.barrierFrame.grid(column=1,row=0,sticky='nsew',padx=5,pady=5)

    self.gaps = [[50,75],[200,250]]
    
    ttk.Button(self.barrierFrame, text='Add Opening').grid(column=0, row=0,sticky='nsew', padx=5, pady=5) #todo: add command

    
#The simulation frame
    self.simFrame = ttk.Labelframe(self.root, text='Simulation')
    self.simFrame.grid(column=0,row=1,columnspan=2,sticky='nsew',padx=5,pady=5)
    
    #buttons to control the simulation
    ttk.Button(self.simFrame, text='Run', command=lambda: self.start()).grid(column=0, row=0,sticky='nsew', padx=5, pady=5) #todo: add command
    ttk.Button(self.simFrame, text='Stop', command=lambda: self.stop()).grid(column=1, row=0,sticky='nsew', padx=5, pady=5) #todo: add command
    ttk.Button(self.simFrame, text='Reset').grid(column=2, row=0,sticky='nsew', padx=5, pady=5) #todo: add command

    #label to show the current frequency
    freqLabel = ttk.Label(self.simFrame, text="Frequency = ")
    freqLabel.grid(column=4, row=1,sticky='nsw', padx=5, pady=5)
    #slider to set the frequency
    ttk.Label(self.simFrame, text="Min").grid(column=3, row=0,sticky='nsw', padx=5, pady=5) #todo: give value for min
    freqScale = ttk.Scale(self.simFrame, orient=Tkinter.HORIZONTAL, from_=self.minFreq, to=self.maxFreq) #todo: units, disable when simulation is running
    #command to update frequency label when slider is changed
    def setFreqLabel(arg=None):
      freqLabel.config(text = "Frequency = " + str(int(freqScale.get())))
    freqScale.config(command = setFreqLabel)
    freqScale.grid(column=4, row=0,sticky='nsew', padx=5, pady=5)
    ttk.Label(self.simFrame, text="Max").grid(column=5, row=0,sticky='nsw', padx=5, pady=5) #todo: give value for min
    #give initial frequency value
    setFreqLabel() #todo: not working; starts at zero
  
  
    self.redrawCanvas();
  
  def updateEzPlot(self):
    #the numpy array; just for testing purposes
    #notes on array
    #1)the only floats Image.fromstring can handle are 32bit, so need to do that conversion
    #2)0 (and below) are black, 255 and above are white, shades of gray inbetween
    #3)need to store array in (height, width) format
    data = np.float32(np.random.rand(200, 400))*255
    im = Image.fromstring('F', (data.shape[1], data.shape[0]), data.tostring())
    self.ezPlot = ImageTk.PhotoImage(image=im) #need to store it so it doesn't get garbage collected, otherwise it won't display correctly on the canvas
  
  def redrawCanvas(self):
    #first, clear everything off the canvas (but don't delete the canvas itself)
    self.canvas.delete('all')
    
    self.updateEzPlot()
    self.canvas.create_image(0,0,image=self.ezPlot,anchor=Tkinter.NW)    
    
    #next, draw the barrier
    #draw the first part
    self.canvas.create_line([(self.barrierX,0),(self.barrierX,self.gaps[0][0])], width=1, fill='red')
    #draw intermediate parts
    if len(self.gaps) > 1:
      for i in range(len(self.gaps)):
	if i < (len(self.gaps)-1):
	  self.canvas.create_line([(self.barrierX,self.gaps[i][1]),(self.barrierX,self.gaps[i+1][0])], width=1, fill='red')    
    #draw the last part
    self.canvas.create_line([(self.barrierX,self.gaps[-1][1]),(self.barrierX,self.viewHeight)], width=1, fill='red')    
    
  def start(self):
    self.cont = True
    self.run()
    
  def stop(self):
    self.cont = False
    
  def run(self):
    if self.cont:
      #todo: step simulation
      self.redrawCanvas()
      self.root.after(1000,self.run) #todo: adjust time  
  
  
  
  
  
  
  
  
  
if __name__ == "__main__":
#the root window
  gui = Interface() #make the interface
  gui.root.mainloop() #set everything in motion