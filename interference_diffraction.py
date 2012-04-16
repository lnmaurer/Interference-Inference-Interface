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
  
  Db = dt/mu0/d
  Cb = dt/epsilon0/d
  
  NEZx = Nx #Ez at x endpoints
  NEZy = Ny #Ez at y endpoints
  #range updated by Yee algorithm:
  EZx_range = slice(barrierX+1, NEZx-1) #right end is RBC and left side is explicitly calculated, so not updated using Yee algorithm
  EZy_range = slice(1, NEZy-1) #top and bottom are RBCs, so not updated using Yee algorithm
  #range updated by explicit formula:
  EZx_ex_range = slice(0, barrierX+1)
  EZy_ex_range = slice(0, NEZy)
  
  NHXx = Nx     #Hx at x endpoints
  NHXy = Ny - 1 #Hx not at y endpoints, so -1
  HXx_range = slice(0, NHXx) #all positions updated using Yee
  HXy_range = slice(0, NHXy)

  NHYx = Nx - 1 #Hy not at x endpoints, so -1
  NHYy = Ny     #Hy at y endpoints
  HYx_range = slice(0, NHYx) #all y positions updated using Yee
  HYy_range = slice(0, NHYy) 
  
  minFreq = 1
  maxFreq = 37

  def __init__(self):
#for the simulation
    self.t = 0 #will hold the current time
    self.cont = False #the simulation isn't currently running
    self.Ez = np.zeros((self.NEZx, self.NEZy,28)) #3rd dimension to keep track of past values of Ez. todo: change it's number to make appropriate to find Ez_rms
    self.Hx = np.zeros((self.NHXx, self.NHXy))
    self.Hy = np.zeros((self.NHYx, self.NHYy))
    
    #for testing only
    self.lamb = 20*self.d
    self.k = 2*np.pi/self.lamb
    self.omega = 2*np.pi*self.c/self.lamb
    self.tau = self.lamb/self.c
    
    print str(self.tau/self.dt)
    
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
    data = np.float32((np.transpose(self.Ez[:,:,0])+1)/2*256)
    im = Image.fromstring('F', (data.shape[1], data.shape[0]), data.tostring())
    self.ezPlot = ImageTk.PhotoImage(image=im) #need to store it so it doesn't get garbage collected, otherwise it won't display correctly on the canvas
    
  def updateEzRMSPlot(self):
    data = np.float32(np.transpose(np.sqrt(np.sum(np.square(self.Ez),axis=2)/28)))*256*(2**0.5) #RMS of Ez
    im = Image.fromstring('F', (data.shape[1], data.shape[0]), data.tostring())
    self.ezRMSPlot = ImageTk.PhotoImage(image=im) #need to store it so it doesn't get garbage collected, otherwise it won't display correctly on the canvas
    
  def redrawCanvas(self):
    #first, clear everything off the canvas (but don't delete the canvas itself)
    self.canvas.delete('all')
    
    self.updateEzPlot()
    self.canvas.create_image(0,0,image=self.ezRMSPlot,anchor=Tkinter.NW)    
    
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
      self.root.after(1,self.run) #todo: adjust time???
  
  def step(self):
    #first, handle the Ez source -- area from the barrier left
    #x and y for points on the barrier and to the left
    x, y = self.d * np.mgrid[self.EZx_ex_range, self.EZy_ex_range]
    tr = self.t - x/self.c
    #want wave to start gradually and propigate at speed of light
    #logistic growth makes it come in gradually and use of retarded time there and in step function at end enforces propigation
    self.Ez[self.EZx_ex_range, self.EZy_ex_range, 0] = np.sin(self.k*x-self.omega*self.t)/(1+np.exp(-(tr-3*self.tau)/self.tau))*((tr > 0).astype(float))

    #now, enforce Ez=0 on barrier
    invGaps = [ypos for pair in self.gaps for ypos in pair] #flatten self.gaps
    invGaps.insert(0,0) #put zero for the first item
    invGaps.append(self.Ny) #put Ny for the last item
    #now invGaps looks like [0, start of first gap, end of first gap, start of second gap, end of second gap,...,Ny]
    for i in range(0,len(invGaps),2):
      self.Ez[self.barrierX, invGaps[i]:invGaps[i+1], 0] = 0     
    
    
    #next, take care of Hx and Hy using the standard Yee algorithm
    self.Hx[self.HXx_range, self.HXy_range] = self.Hx[self.HXx_range, self.HXy_range] + self.Db*(self.Ez[self.HXx_range, self.HXy_range, 0] - self.Ez[self.HXx_range, 1:(self.NHXy+1), 0]) #HXy_range+1
    self.Hy[self.HYx_range, self.HYy_range] = self.Hy[self.HYx_range, self.HYy_range] + self.Db*(self.Ez[1:(self.NHYx+1), self.HYy_range, 0] - self.Ez[self.HYx_range, self.HYy_range, 0]) #HYx_range+1
    
    #move old Ez back in the array
    self.Ez = np.roll(self.Ez,1,axis=2)
    
    #do the normal Yee updates on Ez in the relevant range
    self.Ez[self.EZx_range, self.EZy_range, 0] = self.Ez[self.EZx_range, self.EZy_range, 1] + self.Cb*(self.Hy[self.EZx_range, self.EZy_range] - self.Hy[self.barrierX:(self.NEZx-2),self.EZy_range] + self.Hx[self.EZx_range, 0:self.NEZy-2] - self.Hx[self.EZx_range, self.EZy_range])
    
    
    ##now take care of the Mur RBCs
    Ma = (self.c*self.dt - self.d)/(self.c*self.dt + self.d)
    Mb = 2*self.d/(self.c*self.dt + self.d)
    Mc = (self.c*self.dt)**2/2/self.d/(self.c*self.dt + self.d)
    ##for x=NEZx-1
    rng   = slice(1,self.NEZy-1) #range of everything in y except the corners
    rngp1 = slice(2,self.NEZy)
    rngm1 = slice(0,self.NEZy-2)
    self.Ez[-1,rng,0]  = -self.Ez[-2,rng,2] + Ma*(self.Ez[-2,rng,0] + self.Ez[-1,rng,2]) + Mb*(self.Ez[-1,rng,1] + self.Ez[-2,rng,1]) + Mc*(self.Ez[-1,rngp1,1] - 2*self.Ez[-1,rng,1] + self.Ez[-1,rngm1,1] + self.Ez[-2,rngp1,1] - 2*self.Ez[-2,rng,1] + self.Ez[-2,rngm1,1])

    #for y=0
    rng = slice(self.barrierX+2,self.NEZx-1) #range of everything in x except the corners
    rngp1 = slice(self.barrierX+3,self.NEZx)
    rngm1 = slice(self.barrierX+1,self.NEZx-2)
    self.Ez[rng,0,0]  = -self.Ez[rng,1,2] + Ma*(self.Ez[rng,1,0] + self.Ez[rng,0,2]) + Mb*(self.Ez[rng,0,1] + self.Ez[rng,1,1]) + Mc*(self.Ez[rngp1,0,1] - 2*self.Ez[rng,0,1] + self.Ez[rngm1,0,1] + self.Ez[rngp1,1,1] - 2*self.Ez[rng,1,1] + self.Ez[rngm1,1,1])

    #for y=NEZy
    self.Ez[rng,-1,0]  = -self.Ez[rng,-2,2] + Ma*(self.Ez[rng,-2,0] + self.Ez[rng,-1,2]) + Mb*(self.Ez[rng,-1,1] + self.Ez[rng,-2,1]) + Mc*(self.Ez[rngp1,-1,1] - 2*self.Ez[rng,-1,1] + self.Ez[rngm1,-1,1] + self.Ez[rngp1,-2,1] - 2*self.Ez[rng,-2,1] + self.Ez[rngm1,-2,1])

    #now for the corners
    self.Ez[-1,0,0] = self.Ez[-2,1,2] #bottom right
    self.Ez[-1,-1,0] = self.Ez[-2,-2,2] #top right
    
    ##finially, update the time and the sources todo: don't have to update the sources twice?
    self.t = self.t + self.dt
    self.Ez[self.EZx_ex_range, self.EZy_ex_range, 0] = np.sin(self.k*x-self.omega*self.t)/(1+np.exp(-(tr-3*self.tau)/self.tau))*((tr > 0).astype(float))
  
  
  
  
  
if __name__ == "__main__":
#the root window
  gui = Interface() #make the interface
  gui.root.mainloop() #set everything in motion