# Copyright 2012, 2013, 2023 by the Micromagnum authors.
#
# This file is part of MicroMagnum.
#
# MicroMagnum is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MicroMagnum is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MicroMagnum.  If not, see <http://www.gnu.org/licenses/>.

from magnum.solver import StepHandler
import numpy as np
import math
import os

#Handler for getting the time for plotting
class StorageStepTimeHandler(StepHandler):

    class StorageId(object):
        def __init__(self, var_id, file_fn, field_fn):
            self.__id = var_id
            self.__file_fn = file_fn
            self.__field_fn = field_fn

        def getFileNameForState(self, state):
            return os.path.normpath(self.__file_fn(state))

        def getFieldForState(self, state):
            return self.__field_fn(state)

    def __init__(self, output_dir):
        super(StorageStepTimeHandler, self).__init__()
        self.__ids = {} # maps var_id to StorageId instance
        self.__comments = []
        self.__output_dir = os.path.normpath(output_dir)
        if not os.path.isdir(self.__output_dir): os.makedirs(self.__output_dir)

    def addComment(self, name, fn):
        #if not isinstance(name, str) or not callable(fn): #2to3
        if not isinstance(name, str) or not isinstance(fn, collections.Callable):
            raise TypeError("StorageStepTimeHandler.addComment: 'name' must be a string and 'fn' must be a function")
        self.__comments.append((name, fn))

    def getCommentsForState(self, state):
        return [(name, str(fn(state))) for name, fn in self.__comments]

    def getOutputDirectory(self):
        return self.__output_dir

    def addVariable(self, var_id, file_fn, field_fn = None):
        if not field_fn: field_fn = lambda state: getattr(state, var_id)
        self.__ids[var_id] = StorageStepTimeHandler.StorageId(var_id, file_fn, field_fn)

    def handle(self, state):
        comments = self.getCommentsForState(state)
        for id, sid in self.__ids.items():
            path  = os.path.normpath(self.__output_dir + "/" + sid.getFileNameForState(state))
            field = sid.getFieldForState(state)
            self.store(id, path, field, state.t , comments)

    def store(self, id, path, field, comments):
        raise NotImplementedError("StorageStepTimeHandler.store is purely virtual.")



try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib
    from matplotlib.colors import hsv_to_rgb
    matplotlib.use('Agg')

except:
    class ImageStorage2(StorageStepTimeHandler):
        def __init__(self, output_dir, field_id_or_ids = [], componentIN = "z", componentNameIN="", filetypeIN = "png", plotModeIN = "both", paletteIN = "Red-Green-Blue-Red"):
            raise NotImplementedError("ImageStorage2: Functionality not supported: python 'matplotlib' package could not be loaded.")
else:
    class ImageStorage2(StorageStepTimeHandler):
        def __init__(self, output_dir, field_id_or_ids = [], componentIN = "z", componentNameIN="", filetypeIN = "png", plotModeIN = "both", paletteIN = "Blue-Black-Red", overwrite=False):
            super(ImageStorage2, self).__init__(output_dir)
            component = self.component = componentIN
            componentName = self.componentName = componentNameIN
            filetype = self.filetype = filetypeIN
            palette = self.palette = paletteIN
            plotMode = self.plotMode = plotModeIN
           
            if filetype not in ["png", "jpg", "gif", "eps"]:
                raise ValueError("ImageStorage2: 'filetype' must be either 'png', 'jpg', 'gif' or 'eps'")
            self.__filetype = filetype

            if hasattr(field_id_or_ids, "__iter__") and type(field_id_or_ids) != str:
                field_ids = list(field_id_or_ids)
            else:
                field_ids = [field_id_or_ids]
            if not all(isinstance(x, str) for x in field_ids):
                raise ValueError("ImageStorage2: 'field_id' parameter must be a either a string or a collection of strings.")

            def make_filename_fn(field_id):
                if not overwrite:
                    return lambda state: "%s-%07i.%s" % (field_id, state.step, filetype)
                else:
                    return lambda state: "%s.%s" % (field_id, filetype)


            for field_id in field_ids:
                self.addVariable(field_id, make_filename_fn(field_id))

        def store(self, id, path, field, time, comments):
            writeImage2(self.palette, self.component, self.componentName, self.filetype, field, path, self.plotMode,time=time)

    # quantity in ["xy Angle", "xz Angle", "yz Angle", "x", "y", "z", "abs(x)", "abs(y)", "abs(z)"]:
    # palette in ["Black-White", "White-Black", "Black-White-Black", "White-Black-White", "Blue-Black-Red", "Red-Black-Blue", "Black-Blue", "Blue-Black", "Black-Red", "Red-Black", "Red-Green-Blue-Red"]:
    # filetype in ["eps", "png", "gif", "jpg"]:
    def writeImage2(palette, quantity, quantitystring, filetype, field, figfile, plotMode = "",time=None):
        if plotMode not in ["pixel","both","arrow","colorhsv","colorhsv_arrow"]:
            raise ValueError("ImageStorage2: plotMode must be either 'colorhsv','colorhsv_arrow', 'pixel', 'arrow' or 'both'")
        (nx,ny,nz) = field.mesh.num_nodes
        (dx,dy,dz) = field.mesh.delta
        dx *= 1e9
        dy *= 1e9
        dz *= 1e9

        #keep aspect ratio of simulated area similar to full picture
        pixelsX = 1024
        pixelsY = int(float(ny)/float(nx)*(pixelsX))
        dpi = 150

        if pixelsY < 200: pixelsY = 200
        
        fig,ax = plt.subplots(figsize=(pixelsX/dpi,pixelsY/dpi),dpi=dpi)
        ax.set_xlim(0,nx*dx)
        ax.set_ylim(0,ny*dy)
        ax.set_xlabel(r"$x\ [\mathrm{nm}]$")
        ax.set_ylabel(r'$y\ [\mathrm{nm}]$')
        if time is not None:
           ax.set_title(fr"$t={time*1e9:.4f}\ \mathrm{{ns}}$")
        ax.set_aspect(1)

        if plotMode in ["pixel","both"]:
            cax = make_axes_locatable(ax).append_axes('right','4%',pad='1.5%')

            if palette == "Black-White":
                pal = LinearSegmentedColormap.from_list('',[[0,0,0],[1,1,1]])
            if palette == "White-Black":
                pal = LinearSegmentedColormap.from_list('',[[1,1,1],[0,0,0]])
            if palette == "Black-White-Black":
                pal = LinearSegmentedColormap.from_list('',[[0,0,0],[1,1,1],[0,0,0]])
            if palette == "White-Black-White":
                pal = LinearSegmentedColormap.from_list('',[[1,1,1],[0,0,0],[1,1,1]])
            if palette == "Blue-Black-Red":
                pal = LinearSegmentedColormap.from_list('',[[0,0,1],[0,0,0.5],[0,0,0],[0.5,0,0],[1,0,0]])
            if palette == "Red-Black-Blue":
                pal = LinearSegmentedColormap.from_list('',[[1,0,0],[0.5,0,0],[0,0,0],[0,0,0.5],[0,0,1]])
            if palette == "Black-Blue":
                pal = LinearSegmentedColormap.from_list('',[[0,0,0],[0,0,1]])
            if palette == "Blue-Black":
                pal = LinearSegmentedColormap.from_list('',[[0,0,1],[0,0,0]])
            if palette == "Black-Red":
                pal = LinearSegmentedColormap.from_list('',[[0,0,0],[1,0,0]])
            if palette == "Red-Black":
                pal = LinearSegmentedColormap.from_list('',[[1,0,0],[0,0,0]])
            if palette == "Red-Green-Blue-Red":
                pal = LinearSegmentedColormap.from_list('',[[1,0,0],[0,1,0],[0,0,1],[1,0,0]])
            if palette == "Teal-White-Red":
                pal = LinearSegmentedColormap.from_list('',[[0,0.5,1],[1,1,1],[1,0,0]])
            if palette == "Blue-White-Red": 
                pal = LinearSegmentedColormap.from_list('',[[0,0,1],[0,1,1],[1,1,1],[1,1,0],[1,0,0]])


        if type(field).__name__ == "VectorField":
            #set amount of arrows
            minPixelPerArrow = 30
            pixelsPerCell = float(pixelsX)/float(nx)
            sub = float(minPixelPerArrow)/float(pixelsPerCell)


            sub = max(sub,1)
            #if sub < 2:
            #    subX = round(float(nx)/float(minPixelPerArrow),0)
            #    subY = round(float(ny)/float(minPixelPerArrow),0)
            #    if subX < subY: sub = subX
            #    else: sub = subY
                
            sub=int(sub)


            maxM = 0
            for y in range(0, ny):
                for x in range(0, nx):
                    (Mx, My, Mz) = field.get(x, y, 0)
                    norm = math.sqrt(Mx*Mx+My*My+Mz*Mz)
                    maxM = max(maxM, norm)

            pos = np.zeros((int(ny/sub),int(nx/sub),2))
            vec = np.zeros((int(ny/sub),int(nx/sub),3))

            for y in range(0, int(ny/sub)):
                for x in range(0, int(nx/sub)):
                    sumMx = 0.0
                    sumMy = 0.0
                    sumMz = 0.0
                    num = 0.0
                    for sy in range(0, sub):
                        for sx in range(0,sub):
                            (Mx, My, Mz) = field.get(sub*x+sx, sub*y+sy, 0)
                            sumMx += Mx
                            sumMy += My
                            sumMz += Mz
                            num += 1.0
                    pos[y,x] = np.array([(x+0.5)*dx*sub,(y+0.5)*dy*sub])
                    vec[y,x] = np.array([sumMx/num/maxM,sumMy/num/maxM,sumMz/num/maxM])

            m = np.zeros((ny,nx,3),dtype=np.float32)
            for b in range(0,ny):
                for a in range(0,nx):
                    (Mx,My,Mz) = field.get(a, b, 0)
                    m[b,a] = np.array([Mx/maxM, My/maxM, Mz/maxM])
                    #backgr.append([a, b, Mx/maxM, My/maxM, Mz/maxM])
            


            if plotMode in ["colorhsv","colorhsv_arrow"]:
                #m = np.swapaxes(m,0,1)
                L = ((m[:,:,2]+1)/2)
                H = ((np.arctan2(m[:,:,1],m[:,:,0]))%(2*np.pi)/(2*np.pi))
                V = L+np.min(np.array([L,1-L]),axis=0)
                S = np.zeros_like(L)
                mask = V!=0
                S[mask] = 2*(1-L[mask]/V[mask])
                G = np.array([H,S,V])
                G = np.moveaxis(G,0,2)
                img = hsv_to_rgb(G)
                img[np.linalg.norm(m,axis=2)==0] = np.array([0.5,0.5,0.5])
                ax.imshow(img,extent=(0,nx*dx,0,ny*dy),origin='lower')
            if plotMode == "colorhsv_arrow":
                vecfield = m
                ax.quiver(pos[:,:,0]-vecfield[:,:,0]*0.7*dx*sub/2,pos[:,:,1]-vecfield[:,:,1]*0.7*dy*sub/2,vecfield[:,:,0]*0.7*dx*sub,vecfield[:,:,1]*0.7*dy*sub,scale=1,scale_units='x')  
            elif plotMode == "pixel" or plotMode== "both" or plotMode == "arrow":
                if quantity == "xy Angle":
                    vmin,vmax = -np.pi,np.pi
                    backfield = np.arctan2(m[:,:,0],m[:,:,1])
                    posfield = pos
                    vecfield = vec[:,:,0:2]	
                if quantity == "xz Angle":
                    vmin,vmax = -np.pi,np.pi
                    backfield = np.arctan2(m[:,:,0],m[:,:,2])
                    posfield = pos
                    vecfield = vec[:,:,0:2]		
                if quantity == "yz Angle":
                    vmin,vmax = -np.pi,np.pi
                    backfield = np.arctan2(m[:,:,1],m[:,:,2])
                    posfield = pos
                    vecfield = vec[:,:,0:2]					
                if quantity == "x":
                    vmin,vmax = -1,1
                    backfield = m[:,:,0]
                    posfield = pos
                    vecfield = vec[:,:,0:2]	
                if quantity == "y":
                    vmin,vmax = -1,1
                    backfield = m[:,:,1]
                    posfield = pos
                    vecfield = vec[:,:,0:2]			
                if quantity == "z":
                    vmin,vmax = -1,1
                    backfield = m[:,:,2]
                    posfield = pos
                    vecfield = vec[:,:,0:2]		
                if quantity == "abs(x)":
                    vmin,vmax = 0,1
                    backfield = abs(m[:,:,0])
                    posfield = pos
                    vecfield = vec[:,:,0:2]	
                if quantity == "abs(y)":
                    vmin,vmax = 0,1
                    backfield = abs(m[:,:,1])
                    posfield = pos
                    vecfield = vec[:,:,0:2]	
                if quantity == "abs(z)":
                    vmin,vmax = 0,1
                    backfield = abs(m[:,:,2])
                    posfield = pos
                    vecfield = vec[:,:,0:2]			
                if plotMode == "pixel" or plotMode== "both":
                    c = ax.imshow(backfield,extent=(0,nx*dx,0,ny*dy),vmin=vmin,vmax=vmax,cmap=pal,origin='lower')
                    fig.colorbar(c,cax=cax,label=quantitystring)
                if plotMode == "arrow" or plotMode== "both":
                    ax.quiver(posfield[:,:,0]-vecfield[:,:,0]*0.7*dx*sub/2,posfield[:,:,1]-vecfield[:,:,1]*0.7*dy*sub/2,vecfield[:,:,0]*0.7*dx*sub,vecfield[:,:,1]*0.7*dy*sub,scale=1,scale_units='x')            
        else:
            if plotMode == "pixel" or plotMode == "both":
                m = np.zeros((ny,nx),dtype=np.float32)
                for b in range(0,ny):
                    for a in range(0,nx):
                        m[b,a] = field.get(a, b, 0)

                c = ax.imshow(m,extent=(0,nx*dx,0,ny*dy),cmap=pal,origin='lower')#vmin=vmin,vmax=vmax
                fig.colorbar(c,cax=cax,label=quantitystring)

        fig.tight_layout()
        fig.savefig(figfile)
        plt.close(fig)