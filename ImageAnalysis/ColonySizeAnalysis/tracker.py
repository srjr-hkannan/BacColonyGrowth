class IndexTracker(object):
    def __init__(self, ax, X, plottype, datatype, init_idx=None, args_dict={}):
        self.ax = ax
        self.X = X
        self.plottype = plottype
        self.datatype = datatype
        self.args_dict = args_dict
        
        if init_idx is not None:
            self.ind = init_idx
        else:
            self.ind = self.slices//2

        
        if args_dict.get("title") is not None:
            ax.set_title(args_dict.get("title"))
        else:
            ax.set_title('use scroll wheel to navigate images')
            
        self.zvals = self.args_dict.get("zvals")
        self.zunit = self.args_dict.get("zunit")
        self.zstr = self.args_dict.get("zstr")
        self.zval = None
        if (self.zvals is not None):
            self.zval = self.zvals[self.ind]
        
        if self.plottype == "line_lst":
            self.slices = len(self.X[0])
            xvals = self.X[0][self.ind]
            yvals = self.X[1][self.ind]
            self.im = ax.plot(xvals, yvals, self.datatype)[0]
        else:
            rows, cols, self.slices = X.shape
            if self.plottype == "line":
                # Line plot has data structured as: X[0,:,slice] as the x values, X[1,:,slice] as the y values
                xvals = self.X[0,:,self.ind]
                yvals = self.X[1,:,self.ind]

                self.im = ax.plot(xvals, yvals, self.datatype)[0]
            elif self.plottype == "loglog":
                # Line plot has data structured as: X[0,:,slice] as the x values, X[1,:,slice] as the y values
                xvals = self.X[0,:,self.ind]
                yvals = self.X[1,:,self.ind]

                self.im = ax.loglog(xvals, yvals)[0]

            elif self.plottype == "image":
                # Args dictionary, includes plot parameters
                # x,y limits saved in arg extent

                self.im = ax.imshow(self.X[:, :, self.ind],
                                   extent = self.args_dict.get("extent"))
        self.update()

    def onscroll(self, event):
        print("%s %s" % (event.button, event.step))
        if event.button == "up":
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        if self.plottype == "image":
            self.im.set_data(self.X[:, :, self.ind])
        elif self.plottype == "line_lst":
            xvals = self.X[0][self.ind]
            yvals = self.X[1][self.ind]
            self.im.set_data(xvals,yvals)
        elif self.plottype == "line":
            xvals = self.X[0,:,self.ind]
            yvals = self.X[1,:,self.ind]
            self.im.set_data(xvals, yvals)
        elif self.plottype == "loglog":
            xvals = self.X[0,:,self.ind]
            yvals = self.X[1,:,self.ind]
            self.im.set_data(xvals, yvals)
        # Set zval
        if self.zval is not None:
            self.zval = self.zvals[self.ind]
            self.ax.set_ylabel("%s = %.3f %s" % (self.zstr, self.zval, self.zunit))
        else:
            self.ax.set_ylabel("split %s" % (self.ind))
        self.im.axes.figure.canvas.draw()
