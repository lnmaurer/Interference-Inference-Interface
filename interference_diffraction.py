import Tkinter, tkFileDialog
import ttk
import tkMessageBox
import numpy as np
from PIL import Image, ImageTk

class Interface:
  """The class for the GUI interface"""
  
  mu0      = 1.2566370614e-6    #Vacuum permeability
  epsilon0 = 8.854187817620e-12 # Vacuum permittivity
  c        = 1/(mu0*epsilon0)**0.5
  
  Nx = 800 #width of simulation and view canvas
  Ny = 400 #height of simulation and view canvas
  barrierX = 100 #x position of the barrier
  
  d  = 2.0/Nx
  dt = d/c/2**0.5
  
  NEZx = Nx #Ez at x endpoints
  NEZy = Ny #Ez at y endpoints
  #range updated by Yee algorithm:
  EZx_range = slice(barrierX+2, NEZx-2) #right end is RBC and left side is explicitly calculated, so not updated using Yee algorithm
  EZy_range = slice(1, NEZy-2) #top and bottom are RBCs, so not updated using Yee algorithm
  #range updated by explicit formula:
  EZx_ex_range = slice(0, barrierX+1)
  EZy_ex_range = slice(0, NEZy-1)
  
  NHXx = Nx     #Hx at x endpoints
  NHXy = Ny - 1 #Hx not at y endpoints, so -1
  HXx_range = slice(0, NHXx-1) #all positions updated using Yee
  HXy_range = slice(0, NHXy-1)

  NHYx = Nx - 1 #Hy not at x endpoints, so -1
  NHYy = Ny     #Hy at y endpoints
  HYx_range = slice(0, NHYx-1) #all y positions updated using Yee
  HYy_range = slice(0, NHYy-1) 
  
  minFreq = 1
  maxFreq = 37

  def __init__(self):
#for the simulation
    self.t = 0 #will hold the current time
    self.cont = False #the simulation isn't currently running
    self.Ez = np.zeros((self.Nx, self.Ny))
    
    #for testing only
    self.lamb = 20*self.d
    self.k = 2*np.pi/self.lamb
    self.omega = 2*np.pi*self.c/self.lamb
    self.tau = self.lamb/self.c
    
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
    
    self.canvas = Tkinter.Canvas(self.viewFrame, width=self.Nx, height=self.Ny)
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
    ttk.Button(self.simFrame, text='Run', command=lambda: self.start()).grid(column=0, row=0,sticky='nsew', padx=5, pady=5)
    ttk.Button(self.simFrame, text='Stop', command=lambda: self.stop()).grid(column=1, row=0,sticky='nsew', padx=5, pady=5)
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
    #data = np.float32(np.random.rand(200, 400))*255
    data = np.float32((np.transpose(self.Ez)+1)/2*256)
    im = Image.fromstring('F', (data.shape[1], data.shape[0]), data.tostring())
    self.ezPlot = ImageTk.PhotoImage(image=im) #need to store it so it doesn't get garbage collected, otherwise it won't display correctly on the canvas
  
  def redrawCanvas(self):
    #first, clear everything off the canvas (but don't delete the canvas itself)
    self.canvas.delete('all')
    
    self.updateEzPlot()
    self.canvas.create_image(0,0,image=self.ezPlot,anchor=Tkinter.NW)    
    
    #next, draw the barrier
    invGaps = [ypos for pair in self.gaps for ypos in pair] #flatten self.gaps
    invGaps.insert(0,0) #put zero for the first item
    invGaps.append(self.Ny) #put Ny for the last item
    #now invGaps looks like [0, start of first gap, end of first gap, start of second gap, end of second gap,...,Ny]
    for i in range(0,len(invGaps),2):
      self.canvas.create_line([(self.barrierX,invGaps[i]),(self.barrierX,invGaps[i+1])], width=1, fill='red')    
  
  def reset():
    pass
    #todo: enable frequency slider
    #todo: set t=0
    #todo: clear the canvas
  
  def start(self):
    #todo: reset if it's never been run before
    #todo: don't do anything if it's already running
    #todo: disable frequency slider
    self.cont = True
    self.run()
    
  def stop(self):
    self.cont = False
    
  def run(self):
    if self.cont:
      self.step()
      self.redrawCanvas()
      self.root.after(100,self.run) #todo: adjust time  
  
  def step(self):
    #first, handle area from the barrier left
    #x and y for points on the barrier and to the left
    x, y = self.d * np.mgrid[self.EZx_ex_range, self.EZy_ex_range]
    tr = self.t - x/self.c
    #want wave to start gradually and propigate at speed of light
    #logistic growth makes it come in gradually and use of retarded time there and in step function at end enforces propigation
    self.Ez[self.EZx_ex_range, self.EZy_ex_range] = np.sin(self.k*x-self.omega*self.t)/(1+np.exp(-(tr-3*self.tau)/self.tau))*((tr > 0).astype(float))

    #now, enforce Ez=0 on barrier
    invGaps = [ypos for pair in self.gaps for ypos in pair] #flatten self.gaps
    invGaps.insert(0,0) #put zero for the first item
    invGaps.append(self.Ny) #put Ny for the last item
    #now invGaps looks like [0, start of first gap, end of first gap, start of second gap, end of second gap,...,Ny]
    for i in range(0,len(invGaps),2):
      self.Ez[self.barrierX, invGaps[i]:invGaps[i+1]] = 0     
    
    
    self.t = self.t + self.dt
  
  
  
  
  
if __name__ == "__main__":
#the root window
  gui = Interface() #make the interface
  gui.root.mainloop() #set everything in motion