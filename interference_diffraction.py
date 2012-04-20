import Tkinter, tkFileDialog
import ttk
import tkMessageBox
import time
import math
#import numpy as np
from numpy import *
import numpy
from PIL import Image, ImageTk, ImageDraw, ImageFont

class Interface:
  """The class for the GUI interface"""
  
  mu0      = 1.2566370614e-6    #Vacuum permeability
  epsilon0 = 8.854187817620e-12 # Vacuum permittivity
  c        = 1/(mu0*epsilon0)**0.5
  
  Nx = 600 #width of simulation and view canvas
  Ny = 300 #height of simulation and view canvas
  barrierX = 100 #x position of the barrier
  
  plotD = 100 #dimension of plot
  
  d  = 2.0/Nx #spatial grid element size
  dt = d/c/2**0.5 #time step
  
  tStable = (barrierX + sqrt((Nx-barrierX)**2 + Ny**2))*d/c #time until stability is reached: to barrier then longest way across right domain
  
  #Yee algorithm update coefficients
  Db = dt/mu0/d
  Cb = dt/epsilon0/d
  
  NEZx = Nx #number of Ez grid points in x direction -- goes all the way to the edge
  NEZy = Ny #number of Ez grid points in y direction -- goes all the way to the edge
  #range updated by Yee algorithm:
  EZx_range = slice(barrierX+1, NEZx-1) #right end is RBC and left side is explicitly calculated, so not updated using Yee algorithm
  EZy_range = slice(1, NEZy-1) #top and bottom are RBCs, so not updated using Yee algorithm
  #range updated by explicit formula:
  EZx_ex_range = slice(0, barrierX+1)
  EZy_ex_range = slice(0, NEZy)
  #everything
  EZx_all = slice(0,NEZx)
  EZy_all = slice(0,NEZy)
  
  NHXx = Nx     #Hx goes all the way to x endpoints
  NHXy = Ny - 1 #Hx not at y endpoints, so -1
  HXx_range = slice(0, NHXx) #all positions updated using Yee
  HXy_range = slice(0, NHXy)

  NHYx = Nx - 1 #Hy not at x endpoints, so -1
  NHYy = Ny     #Hy goes all the way to y endpoints
  HYx_range = slice(0, NHYx) #all y positions updated using Yee
  HYy_range = slice(0, NHYy) 

  def __init__(self):
#for the simulation
    self.t = 0 #will hold the current time
    self.tAveraging = 0 #hold the time since averaging started
    self.running = False #the simulation isn't currently running
    self.Ez = zeros((self.NEZx, self.NEZy,3)) #3rd dimension to keep track of past values of Ez.
    self.EzSQ = zeros((self.NEZx, self.NEZy))
    self.Hx = zeros((self.NHXx, self.NHXy))
    self.Hy = zeros((self.NHYx, self.NHYy))
    self.maxY = 1 #contains the largest Ez seen
    
    #for testing only???
    self.lamb = 20*self.d
    self.k = 2*pi/self.lamb
    self.omega = 2*pi*self.c/self.lamb
    self.tau = self.lamb/self.c
#misc
    self.font = ImageFont.truetype("BebasNeue.otf", 90)
    self.haveRestartedAvg = False
    self.fastForwarding = False
#The root window
    self.root = Tkinter.Tk()
    self.root.title("Leon's New Fangled Interference & Diffraction Simulator")
    
#The menubar and menus
    menubar = Tkinter.Menu(self.root)
  
    #the file menu
    filemenu = Tkinter.Menu(menubar, tearoff=0)
    #filemenu.add_command(label="Save Experiment As", accelerator="Ctrl+S")
    #filemenu.add_command(label="Load Experiment", accelerator="Ctrl+O")
    #filemenu.add_separator()
    filemenu.add_command(label="Exit", accelerator="Ctrl+Q", command=self.root.quit)
    #bind keys to the actions
    self.root.bind_all('<Control-q>', lambda arg: self.root.quit())
    menubar.add_cascade(label="File", menu=filemenu)
  
    #the edit menu
    editmenu = Tkinter.Menu(menubar, tearoff=0)
    editmenu.add_command(label="Cut", accelerator="Ctrl+X") #todo: add command
    editmenu.add_command(label="Copy", accelerator="Ctrl+C") #todo: add command
    editmenu.add_command(label="Paste", accelerator="Ctrl+V") #todo: add command
    menubar.add_cascade(label="Edit", menu=editmenu)
    
    #the view menu
    viewmenu = Tkinter.Menu(menubar, tearoff=0)
    self.traceSetting = Tkinter.StringVar(value="line")
    viewmenu.add_radiobutton(label="Trace Setting:", state="disabled")
    viewmenu.add_radiobutton(label="Lines", variable=self.traceSetting, value="line", command=self.conditionalRedraw)
    viewmenu.add_radiobutton(label="Dots", variable=self.traceSetting, value="dot", command=self.conditionalRedraw)
    
    viewmenu.add_separator()
    self.avgSetting = Tkinter.StringVar(value='amp')
    viewmenu.add_radiobutton(label="Displayed Average:", state="disabled")
    viewmenu.add_radiobutton(label="Amplitude", variable=self.avgSetting, value="amp", command=self.conditionalRedraw)
    viewmenu.add_radiobutton(label="Ez_rms", variable=self.avgSetting, value="rms", command=self.conditionalRedraw)
    viewmenu.add_radiobutton(label="Ez_rms^2", variable=self.avgSetting, value="sq", command=self.conditionalRedraw)
    menubar.add_cascade(label="View", menu=viewmenu)
    
    #the simulation menu
    simmenu = Tkinter.Menu(menubar, tearoff=0)
    simmenu.add_command(label="Run", accelerator="Ctrl+R", command=self.start)
    self.root.bind("<Control-r>", lambda arg: self.start())
    simmenu.add_command(label="Stop", accelerator="Ctrl+S", command=self.stop)
    self.root.bind("<Control-s>", lambda arg: self.stop())
    simmenu.add_command(label="Reset Simulation", accelerator="Ctrl+Shift+R", command=self.reset)
    self.root.bind("<Control-R>", lambda arg: self.reset())
    simmenu.add_command(label="Reset Average", accelerator="Ctrl+Shift+A", command=self.resetIntensity)
    self.root.bind("<Shift-A>", lambda arg: self.resetIntensity())
    simmenu.add_command(label="Fast Forward", accelerator="Ctrl+F", command=self.fastForward)
    self.root.bind("<Control-f>", lambda arg: self.fastForward())
    menubar.add_cascade(label="Simulation", menu=simmenu)
    
    self.root.config(menu=menubar)
    
#The view frame
    self.viewFrame = ttk.Labelframe(self.root, text='View')
    self.viewFrame.grid(column=0, row=0, sticky='nsew',padx=5,pady=5)
    
    self.Ezcanvas = Tkinter.Canvas(self.viewFrame, width=self.Nx, height=self.Ny)
    self.Ezcanvas.grid(column=0, row=0, columnspan=3, rowspan=3, sticky='nsew', padx=5, pady=5)    

    self.HorizPlotCanvas = Tkinter.Canvas(self.viewFrame, width=self.Nx, height=self.plotD)
    self.HorizPlotCanvas.grid(column=0, row=3, columnspan=3, rowspan=5, sticky='nsew', padx=5, pady=5)
    
    self.xStringVar = Tkinter.StringVar()
    self.yStringVar = Tkinter.StringVar()
    self.tStringVar = Tkinter.StringVar()
    self.EzStringVar = Tkinter.StringVar()
    self.EzRMSStringVar = Tkinter.StringVar()
    ttk.Label(self.viewFrame, textvariable=self.xStringVar).grid(column=3, row=3, sticky='nsew', padx=5, pady=0)
    ttk.Label(self.viewFrame, textvariable=self.yStringVar).grid(column=3, row=4, sticky='nsew', padx=5, pady=0)
    ttk.Label(self.viewFrame, textvariable=self.tStringVar).grid(column=3, row=5, sticky='nsew', padx=5, pady=0)    
    ttk.Label(self.viewFrame, textvariable=self.EzStringVar).grid(column=3, row=6, sticky='nsew', padx=5, pady=0)    
    ttk.Label(self.viewFrame, textvariable=self.EzRMSStringVar).grid(column=3, row=7, sticky='nsew', padx=5, pady=0)    
    
    self.VertPlotCanvas1 = Tkinter.Canvas(self.viewFrame, width=self.plotD, height=self.Ny)
    self.VertPlotCanvas1.grid(column=3, row=0, columnspan=1, rowspan=3, sticky='nsew', padx=5, pady=5)
    self.VertPlotCanvas2 = Tkinter.Canvas(self.viewFrame, width=self.plotD, height=self.Ny)
    self.VertPlotCanvas2.grid(column=3, row=8, columnspan=1, rowspan=3, sticky='nsew', padx=5, pady=5)
    
    self.sliceY = 60 #position of horizontal slice todo: put in middle, set to 60 for testing
    self.sliceX = self.Nx/2
    
    self.EzRMScanvas = Tkinter.Canvas(self.viewFrame, width=self.Nx, height=self.Ny)
    self.EzRMScanvas.grid(column=0, row=8, columnspan=3, rowspan=3, sticky='nsew', padx=5, pady=5)
   
    #make it so that, after dragging an element has ceased, the binding is reset so that further dragging won't move the element unless it gets clicked again first
    #we do this by binding to any motion on any canvas
    self.Ezcanvas.bind("<Motion>", self.clearCanvasBindings)
    self.HorizPlotCanvas.bind("<Motion>", self.clearCanvasBindings)
    self.EzRMScanvas.bind("<Motion>", self.clearCanvasBindings)
    
    #make the arrow keys control the slice positions
    self.root.bind("<Up>", lambda arg: self.setSliceY(self.sliceY-1))
    self.root.bind("<Down>", lambda arg: self.setSliceY(self.sliceY+1))
    self.root.bind("<Right>", lambda arg: self.setSliceX(self.sliceX+1))
    self.root.bind("<Left>", lambda arg: self.setSliceX(self.sliceX-1))
    
#The barrier frame frame and intial barrier setup
    self.barrierFrame = ttk.Labelframe(self.root, text='Barrier')
    self.barrierFrame.grid(column=1,row=0,sticky='nsew',padx=5,pady=5)

    ttk.Button(self.barrierFrame, text='Add Opening', command=self.addOpening).grid(column=0, row=0, columnspan=2, sticky='nsew', padx=5, pady=5)
    ttk.Label(self.barrierFrame, text="Bottom:").grid(column=0, row=1, sticky='nes', padx=5, pady=5)
    self.bottomEntry = Tkinter.Spinbox(self.barrierFrame, width=4, from_=0, to=self.Nx)
    self.bottomEntry.grid(column=1, row=1, sticky='nsw', padx=5, pady=5)
    ttk.Label(self.barrierFrame, text="Top:").grid(column=0, row=2, sticky='nes', padx=5, pady=5)
    self.topEntry = Tkinter.Spinbox(self.barrierFrame, width=4, from_=0, to=self.Nx)
    self.topEntry.grid(column=1, row=2, sticky='nsw', padx=5, pady=5)
    
    self.gaps = [[50,70],[230,250]]
    self.barrierFrames = []
    
    self.redrawBarrierFrame()

#almost done
    self.redrawCanvases();    
  
  def addOpening(self):
    bot = int(self.bottomEntry.get())
    top = int(self.topEntry.get())
    if bot > 0 and top < self.Ny and bot < top: #some initial checks that top and bot better pass
      newOrder = self.gaps[:]
      newOrder.append([bot,top])
      newOrder.sort()
      newOrderFlat = [ypos for pair in newOrder for ypos in pair] #flatten
      if (newOrderFlat == sorted(newOrderFlat) #catches things like [50,100],[90,120] -- overlapping openings
         and len(list(set(newOrderFlat))) == len(newOrderFlat)): #catchs things like [50,100],[100,120] -- openings with no space between them
	self.gaps = newOrder
	self.redrawBarrierFrame()
	self.conditionalRedraw()
	
  
  def redrawBarrierFrame(self):
    self.gaps.sort()
    for b in self.barrierFrames:
      b.destroy()
    self.barrierFrames = []
    self.strVars = [] #need to save StringVars or else they get garbage collected
    self.distances = []
    r = 3 #rows 0,1,2 already taken by widgets for adding an opeing
    
    for gap in self.gaps:
      bn = r-3 #barrier number
      frame = ttk.Labelframe(self.barrierFrame, text="Opening {}".format(bn))
      frame.grid(column=0, row=r, columnspan=2, sticky='nesw', padx=5, pady=5)
      
      top = Tkinter.IntVar()
      bottom = Tkinter.IntVar()
      distance = Tkinter.StringVar()
      self.setDistance(gap, distance)
      ttk.Label(frame, text="Top:").grid(column=0, row=0, sticky='nes', padx=5, pady=5)
      #having the following work is kind of tricky; the default parameter in the lambda is critical. See <http://mail.python.org/pipermail/tutor/2005-November/043360.html>
      entry = Tkinter.Spinbox(frame, width=4, textvariable=top, from_=0, to=self.Nx, command=lambda n=bn, tv=top: self.updateBarrierTop(n,tv))
      entry.bind("<Return>",lambda arg, n=bn, tv=top: self.updateBarrierTop(n,tv))
      entry.grid(column=1, row=0, sticky='nsw', padx=5, pady=5)
      ttk.Label(frame, text="Bottom:").grid(column=0, row=1, sticky='nes', padx=5, pady=5)
      entry = Tkinter.Spinbox(frame, width=4, textvariable=bottom, from_=0, to=self.Nx, command=lambda n=bn, tv=bottom: self.updateBarrierBottom(n,tv))
      entry.bind("<Return>",lambda arg, n=bn, tv=bottom: self.updateBarrierBottom(n,tv))
      entry.grid(column=1, row=1, sticky='nsw', padx=5, pady=5)
      ttk.Label(frame, textvariable=distance).grid(column=0, row=2, sticky='nes', padx=5, pady=5)
      ttk.Button(frame, text='Remove', command=lambda n=bn: self.removeBarrier(n)).grid(column=0, row=3, sticky='nsew', columnspan=2, padx=5, pady=5)
      top.set(gap[1])
      bottom.set(gap[0])
      self.strVars.append(top)
      self.strVars.append(bottom) 
      self.distances.append(distance)
      
      self.barrierFrames.append(frame)
      r = r + 1
      
  def updateBarrierTop(self, barrierNumber, intVar):
    value = intVar.get()
    if ((barrierNumber == (len(self.gaps)-1)) and (value < self.Ny) and (value > self.gaps[-1][0])) or ((barrierNumber < (len(self.gaps)-1)) and (value < self.gaps[barrierNumber+1][0]) and (value > self.gaps[barrierNumber][0])):
      self.gaps[barrierNumber][1] = value
      self.conditionalRedraw()
    else:
      intVar.set(self.gaps[barrierNumber][1])

  def updateBarrierBottom(self, barrierNumber, intVar):
    value = intVar.get()
    if ((barrierNumber == 0) and (value > 0) and (value < self.gaps[0][1])) or ((barrierNumber > 0) and (value > self.gaps[barrierNumber-1][1]) and (value < self.gaps[barrierNumber][1])):
      self.gaps[barrierNumber][0] = value
      self.conditionalRedraw()
    else:
      intVar.set(self.gaps[barrierNumber][0])    
  
  def setDistance(self, gap, textVar):
    textVar.set(str(int(round(sqrt(((gap[0]+gap[1])/2.0-self.sliceY)**2+(self.barrierX-self.sliceX)**2)))))
  
  def removeBarrier(self, barrierNumber):
    del self.gaps[barrierNumber]
    self.redrawBarrierFrame()
    self.conditionalRedraw()
    
  def clearCanvasBindings(self, eventObj):
    self.Ezcanvas.bind("<B1-Motion>", lambda e: None)
    self.HorizPlotCanvas.bind("<B1-Motion>", lambda e: None)
    self.EzRMScanvas.bind("<B1-Motion>", lambda e: None)
  
  def updateEzPlot(self):
    #the numpy array; just for testing purposes
    #notes on array
    #1)the only floats Image.fromstring can handle are 32bit, so need to do that conversion
    #2)0 (and below) are black, 255 and above are white, shades of gray inbetween
    #3)need to store array in (height, width) format
    data = float32((transpose(self.Ez[:,:,0])/self.maxY+1)/2*256)
    im = Image.fromstring('F', (data.shape[1], data.shape[0]), data.tostring())
    self.ezPlot = ImageTk.PhotoImage(image=im) #need to store it so it doesn't get garbage collected, otherwise it won't display correctly on the canvas
    
  def updateEzRMSPlot(self):
    if self.avgSetting.get() == 'sq':
      data = 256*(float32(transpose(self.EzSQ)/self.tAveraging/self.maxRMSY**2))  
    else:
      data = 256*(float32(transpose(self.EzRMS)/self.maxRMSY))    
    im = Image.fromstring('F', (data.shape[1], data.shape[0]), data.tostring())
    if self.t < self.tStable and self.running:
      count = int(math.ceil((self.tStable - self.t)/self.dt))
      draw = ImageDraw.Draw(im)
      draw.text((self.Nx/3, self.Ny/3), str(count), font=self.font, fill=255.0)
    self.ezRMSPlot = ImageTk.PhotoImage(image=im) #need to store it so it doesn't get garbage collected, otherwise it won't display correctly on the canvas

  def applyAveraging(self, xS, yS, invert=False, othercoord=None):
    if self.avgSetting.get() == 'amp':
      v = int_(numpy.round((1+self.EzRMS[xS,yS]*sqrt(2)/self.maxY)*50))
      if othercoord != None: #in this case, let's plot +/-amplitude
	othercoord.extend(othercoord[::-1]) #extend the coordinates so they're like [0,1...,Nx-1,Nx,Nx,Nx-1,...,1,0]
	v = concatenate((v,100-v[::-1])) #extend the values so they're like [V0,...,Vn,-Vn,...,-V0]
    elif self.avgSetting.get() == 'rms':
      v = int_(numpy.round((self.EzRMS[xS,yS]/self.maxRMSY)*99))
    else:
      v = int_(numpy.round((self.EzSQ[xS,yS]/self.tAveraging/self.maxRMSY**2)*99))
    
    if invert:
      v = 100 - v
    return v
    
  def plot(self, draw, x, y, color):
    if self.traceSetting.get() == 'line':
      draw.line(zip(x,y), fill=color)
    else:
      draw.point(zip(x,y), fill=color)
      
  def updateHorizPlot(self):
    im = Image.new('RGB', (self.Nx,self.plotD))
    draw = ImageDraw.Draw(im)
    x = range(0,self.Nx)
    
    #plot EzRMS
    y = self.applyAveraging(self.EZx_all, self.sliceY, invert=True, othercoord=x)
    self.plot(draw, x, y, 'green')
      
    #plot Ez
    y = int_(numpy.round((-self.Ez[:,self.sliceY,0]/self.maxY+1)*50))
    self.plot(draw, x, y, 'yellow')
      
    self.horizPlot = ImageTk.PhotoImage(image=im)

  def updateVertPlot(self):
    im = Image.new('RGB', (self.plotD, self.Ny))
    draw = ImageDraw.Draw(im)
    y = range(0,self.Ny)
    
    #plot EzRMS
    x = self.applyAveraging(self.sliceX, self.EZy_all, othercoord=y)
    self.plot(draw, x, y, 'green')
      
    #plot Ez
    x = int_(numpy.round((self.Ez[self.sliceX,:,0]/self.maxY+1)*50))
    self.plot(draw, x, y, 'yellow')
      
    self.vertPlot = ImageTk.PhotoImage(image=im)  
  
    
  def redrawCanvases(self):
    #first, clear everything off the canvases (but don't delete the canvases themselves)
    self.Ezcanvas.delete('all')
    self.HorizPlotCanvas.delete('all')
    self.EzRMScanvas.delete('all')
    self.VertPlotCanvas1.delete('all')
    self.VertPlotCanvas2.delete('all')
    
    self.EzRMS = sqrt(self.EzSQ/(self.tAveraging+self.dt/1e6))
    self.maxRMSY = numpy.max(self.EzRMS)
    if self.maxRMSY == 0:
      self.maxRMSY = 1
      
    tempMaxY = numpy.max(abs(self.Ez))
    if tempMaxY > self.maxY:
      self.maxY = tempMaxY
    #print self.maxY
    
    self.updateEzPlot()
    self.Ezcanvas.create_image(0,0,image=self.ezPlot,anchor=Tkinter.NW)
    self.updateHorizPlot()
    self.HorizPlotCanvas.create_image(0,0,image=self.horizPlot,anchor=Tkinter.NW)  
    self.updateEzRMSPlot()
    self.EzRMScanvas.create_image(0,0,image=self.ezRMSPlot,anchor=Tkinter.NW)
    self.HorizPlotCanvas.create_image(0,0,image=self.horizPlot,anchor=Tkinter.NW)  
    self.updateVertPlot()
    self.VertPlotCanvas1.create_image(0,0,image=self.vertPlot,anchor=Tkinter.NW)  
    self.VertPlotCanvas2.create_image(0,0,image=self.vertPlot,anchor=Tkinter.NW)  
    
    #next, draw the barrier
    invGaps = [ypos for pair in self.gaps for ypos in pair] #flatten self.gaps
    invGaps.insert(0,0) #put zero for the first item
    invGaps.append(self.Ny) #put Ny for the last item
    #now invGaps looks like [0, start of first gap, end of first gap, start of second gap, end of second gap,...,Ny]
    for i in range(0,len(invGaps),2):
      self.Ezcanvas.create_line([(self.barrierX,invGaps[i]),(self.barrierX,invGaps[i+1])], width=1, fill='red')    
      self.EzRMScanvas.create_line([(self.barrierX,invGaps[i]),(self.barrierX,invGaps[i+1])], width=1, fill='red')
      
    #now, draw the horizontal slice
    lineID = self.Ezcanvas.create_line([(0,self.sliceY),(self.Nx,self.sliceY)], width=1, fill='yellow', dash='-') 
    self.Ezcanvas.tag_bind(lineID, "<Button-1>",  self.horizClickMethod)
    lineID = self.EzRMScanvas.create_line([(0,self.sliceY),(self.Nx,self.sliceY)], width=1, fill='green', dash='-')
    self.EzRMScanvas.tag_bind(lineID, "<Button-1>",  self.horizClickMethod)    
    lineID = self.VertPlotCanvas1.create_line([(0,self.sliceY),(self.plotD,self.sliceY)], width=1, fill='blue', dash='-')
    self.VertPlotCanvas1.tag_bind(lineID, "<Button-1>",  self.horizClickMethod)
    lineID = self.VertPlotCanvas2.create_line([(0,self.sliceY),(self.plotD,self.sliceY)], width=1, fill='blue', dash='-')
    self.VertPlotCanvas2.tag_bind(lineID, "<Button-1>",  self.horizClickMethod)    
    
    #the vertical slice
    lineID = self.Ezcanvas.create_line([(self.sliceX,0),(self.sliceX,self.Ny)], width=1, fill='yellow', dash='-') 
    self.Ezcanvas.tag_bind(lineID, "<Button-1>",  self.vertClickMethod)
    lineID = self.EzRMScanvas.create_line([(self.sliceX,0),(self.sliceX,self.Ny)], width=1, fill='green', dash='-')
    self.EzRMScanvas.tag_bind(lineID, "<Button-1>",  self.vertClickMethod)
    lineID = self.HorizPlotCanvas.create_line([(self.sliceX,0),(self.sliceX,self.plotD)], width=1, fill='blue', dash='-')
    self.HorizPlotCanvas.tag_bind(lineID, "<Button-1>",  self.vertClickMethod)
    
    #finially, show the corridnates of the slices
    self.xStringVar.set("x=" + str(self.sliceX) + "d")
    self.yStringVar.set("y=" + str(self.sliceY) + "d")
    self.tStringVar.set("t=" + str(int(self.t/self.dt)) + "dt")
    self.EzStringVar.set("Ez={:+.4f}".format(self.Ez[self.sliceX,self.sliceY,0]))
    self.EzRMSStringVar.set("EzRMS={:.4f}".format(self.EzRMS[self.sliceX,self.sliceY]))
    
  def horizClickMethod(self, eventObj):
    self.Ezcanvas.bind('<B1-Motion>', self.horizDragMethod)
    self.EzRMScanvas.bind('<B1-Motion>', self.horizDragMethod)
    self.VertPlotCanvas1.bind('<B1-Motion>', self.horizDragMethod)
    self.VertPlotCanvas2.bind('<B1-Motion>', self.horizDragMethod)
    
  def horizDragMethod(self, eventObj):
    self.setSliceY(eventObj.y)

  def conditionalRedraw(self):
    if not self.running:
      self.redrawCanvases()
      
  def setSliceY(self, y):
    if (y >= 0) and (y < self.Ny):
      self.sliceY = y
      for gap, dist in zip(self.gaps, self.distances):
	self.setDistance(gap, dist)
      self.conditionalRedraw()
    
  def vertClickMethod(self, eventObj):
    self.Ezcanvas.bind('<B1-Motion>', self.vertDragMethod)
    self.EzRMScanvas.bind('<B1-Motion>', self.vertDragMethod)
    self.HorizPlotCanvas.bind('<B1-Motion>', self.vertDragMethod)
    
  def vertDragMethod(self, eventObj):
    self.setSliceX(eventObj.x)
    
  def setSliceX(self, x):
    if (x >= 0) and (x < self.Nx):
      self.sliceX = x
      for gap, dist in zip(self.gaps, self.distances):
	self.setDistance(gap, dist)
      self.conditionalRedraw()  
    
  def resetIntensity(self):
    self.tAveraging = 0
    self.EzSQ = zeros((self.NEZx, self.NEZy))
      
  def reset(self):
    self.t = 0
    self.haveRestartedAvg = False
    self.resetIntensity()
    self.Ez = zeros((self.NEZx, self.NEZy, 3))
    self.maxY = 1.0
    self.conditionalRedraw()
  
  def start(self):
    #todo: reset if it's never been run before
    #todo: don't do anything if it's already running
    #todo: disable frequency slider
    self.running = True
    self.run()
    
  def stop(self):
    self.running = False
  
  def fastForward(self):
    self.tEnd = self.t + self.tStable
    self.fastForwarding = True
    self.root.after(1,self.fastForwardStep)
    
  def fastForwardStep(self):
    if self.t < self.tEnd:
      self.Ezcanvas.delete('all')
      self.EzRMScanvas.delete('all')
      im = Image.new('RGB', (self.Nx,self.Ny))
      draw = ImageDraw.Draw(im)
      draw.text((self.Nx/3, self.Ny/3), str(int((self.tEnd-self.t)/self.dt)), font=self.font, fill='red')
      self.ezPlot = ImageTk.PhotoImage(image=im)
      self.Ezcanvas.create_image(0,0,image=self.ezPlot,anchor=Tkinter.NW)
      self.EzRMScanvas.create_image(0,0,image=self.ezPlot,anchor=Tkinter.NW)
      self.step(avg=False)
      self.root.after(1,self.fastForwardStep)
    else:
      self.fastForwarding = False
      self.resetIntensity()
      self.start()
  
  def run(self):
    if self.running:
      t = time.clock()
      if self.t > self.tStable and not self.haveRestartedAvg: #todo: clean up reset
	self.resetIntensity()
	self.haveRestartedAvg = True
      self.step()
      self.redrawCanvases()
      print str(time.clock()-t)
      if not self.fastForwarding:
	self.root.after(1,self.run)
  
  def step(self, avg=True):
    #so that we don't need a bazillion "self."s
    Ez = self.Ez
    Hx = self.Hx
    Hy = self.Hy
    
    #first, handle the Ez source -- area from the barrier left
    #x and y for points on the barrier and to the left
    x, y = self.d * mgrid[self.EZx_ex_range, self.EZy_ex_range]
    tr = self.t - x/self.c
    #want wave to start gradually and propigate at speed of light
    #logistic growth makes it come in gradually and use of retarded time there and in step function at end enforces propigation
    Ez[self.EZx_ex_range, self.EZy_ex_range, 0] = sin(self.k*x-self.omega*self.t)/(1+exp(-(tr-3*self.tau)/self.tau))*((tr > 0).astype(float))

    #now, enforce Ez=0 on barrier
    invGaps = [ypos for pair in self.gaps for ypos in pair] #flatten self.gaps
    invGaps.insert(0,0) #put zero for the first item
    invGaps.append(self.Ny) #put Ny for the last item
    #now invGaps looks like [0, start of first gap, end of first gap, start of second gap, end of second gap,...,Ny]
    for i in range(0,len(invGaps),2):
      Ez[self.barrierX, invGaps[i]:invGaps[i+1], 0] = 0     
    
    
    #next, take care of Hx and Hy using the standard Yee algorithm
    Hx[self.HXx_range, self.HXy_range] = Hx[self.HXx_range, self.HXy_range] + self.Db*(Ez[self.HXx_range, self.HXy_range, 0] - Ez[self.HXx_range, 1:(self.NHXy+1), 0]) #HXy_range+1
    Hy[self.HYx_range, self.HYy_range] = Hy[self.HYx_range, self.HYy_range] + self.Db*(Ez[1:(self.NHYx+1), self.HYy_range, 0] - Ez[self.HYx_range, self.HYy_range, 0]) #HYx_range+1
    
    #move old Ez back in the array
    Ez = roll(Ez,1,axis=2)
    
    #do the normal Yee updates on Ez in the relevant range
    Ez[self.EZx_range, self.EZy_range, 0] = Ez[self.EZx_range, self.EZy_range, 1] + self.Cb*(Hy[self.EZx_range, self.EZy_range] - Hy[self.barrierX:(self.NEZx-2),self.EZy_range] + Hx[self.EZx_range, 0:self.NEZy-2] - Hx[self.EZx_range, self.EZy_range])
    
    
    ##now take care of the Mur RBCs
    Ma = (self.c*self.dt - self.d)/(self.c*self.dt + self.d)
    Mb = 2*self.d/(self.c*self.dt + self.d)
    Mc = (self.c*self.dt)**2/2/self.d/(self.c*self.dt + self.d)
    ##for x=NEZx-1
    rng   = slice(1,self.NEZy-1) #range of everything in y except the corners
    rngp1 = slice(2,self.NEZy)
    rngm1 = slice(0,self.NEZy-2)
    Ez[-1,rng,0]  = -Ez[-2,rng,2] + Ma*(Ez[-2,rng,0] + Ez[-1,rng,2]) + Mb*(Ez[-1,rng,1] + Ez[-2,rng,1]) + Mc*(Ez[-1,rngp1,1] - 2*Ez[-1,rng,1] + Ez[-1,rngm1,1] + Ez[-2,rngp1,1] - 2*Ez[-2,rng,1] + Ez[-2,rngm1,1])

    #for y=0
    rng = slice(self.barrierX+2,self.NEZx-1) #range of everything in x except right corner
    rngp1 = slice(self.barrierX+3,self.NEZx)
    rngm1 = slice(self.barrierX+1,self.NEZx-2)
    Ez[rng,0,0]  = -Ez[rng,1,2] + Ma*(Ez[rng,1,0] + Ez[rng,0,2]) + Mb*(Ez[rng,0,1] + Ez[rng,1,1]) + Mc*(Ez[rngp1,0,1] - 2*Ez[rng,0,1] + Ez[rngm1,0,1] + Ez[rngp1,1,1] - 2*Ez[rng,1,1] + Ez[rngm1,1,1])

    #for y=NEZy
    Ez[rng,-1,0]  = -Ez[rng,-2,2] + Ma*(Ez[rng,-2,0] + Ez[rng,-1,2]) + Mb*(Ez[rng,-1,1] + Ez[rng,-2,1]) + Mc*(Ez[rngp1,-1,1] - 2*Ez[rng,-1,1] + Ez[rngm1,-1,1] + Ez[rngp1,-2,1] - 2*Ez[rng,-2,1] + Ez[rngm1,-2,1])

    #now for the corners
    Ez[-1,0,0] = Ez[-2,1,2] #bottom right
    Ez[-1,-1,0] = Ez[-2,-2,2] #top right
    
    ##finially, update the time and the sources todo: don't have to update the sources twice?
    self.t = self.t + self.dt
    self.tAveraging = self.tAveraging + self.dt
    Ez[self.EZx_ex_range, self.EZy_ex_range, 0] = sin(self.k*x-self.omega*self.t)/(1+exp(-(tr-3*self.tau)/self.tau))*((tr > 0).astype(float))
    invGaps = [ypos for pair in self.gaps for ypos in pair] #flatten self.gaps
    invGaps.insert(0,0) #put zero for the first item
    invGaps.append(self.Ny) #put Ny for the last item
    #now invGaps looks like [0, start of first gap, end of first gap, start of second gap, end of second gap,...,Ny]
    for i in range(0,len(invGaps),2):
      Ez[self.barrierX, invGaps[i]:invGaps[i+1], 0] = 0  
    
    #strictly speaking, the following will blow up to infinity if you integrate forever (since the FT of a sinusoid is a delta function)
    #however, we're not that patient. plus, floating point limitations will prevent it (once the numbers are large enough, adding a small number to them won't change them)
    if avg:
      self.EzSQ = self.EzSQ + square(Ez[:,:,0])*self.dt
    self.Ez = Ez
    self.Hx = Hx
    self.Hy = Hy  
  
  
  
  
if __name__ == "__main__":
#the root window
  gui = Interface() #make the interface
  gui.root.mainloop() #set everything in motion