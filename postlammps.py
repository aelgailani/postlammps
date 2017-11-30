import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import griddata

class XYZ:
		
	"""docstring for ClassName"""
	def __init__(self, f, 
		columns = ['id', 'type', 'x', 'y', 'radius', 'fx', 'fy', 'Sx', 'Sy', 'Sxy', 'Sz', 'Sxz', 'Syz']):
		
		self.index_rows=[]
		self.headlines_rows=[]
		self.columns_names = columns
		lookup = 'ITEM: ATOMS id type'
		self.content = open(f, 'r').readlines()

		for num, line in enumerate(self.content):
			if lookup in line:
				self.index_rows= self.index_rows + [num]
				self.headlines_rows=  self.headlines_rows + [int(i) for i in range(num-8,num+1)] 

		self.N = [int(s) for s in self.content[3].split() if s.isdigit()]
		
		self.N = int(self.N[0])

		self.df_all = pd.read_csv( f,
                                  skiprows=self.headlines_rows, sep='\s', names=self.columns_names,
                                  index_col='id',engine='python')

		
		self.frames= [ii for ii in range(len(self.index_rows))]
		
		self.df = []
		self.xlo = []
		self.xhi = []
		self.ylo = []
		self.yhi = []
		self.Lx = []
		self.Ly = []
		for ii in self.frames:
			self.df= self.df + [self.df_all.iloc[ii*self.N:ii*self.N+self.N]]

			self.xlo = self.xlo + [self.df[ii]['x'].min(axis=0)]
			self.xhi = self.xhi + [self.df[ii]['x'].max(axis=0)]
			self.ylo = self.ylo + [self.df[ii]['y'].min(axis=0)]
			self.yhi = self.yhi + [self.df[ii]['y'].max(axis=0)]

			self.Lx = self.Lx + [ self.xhi[ii]- self.xlo[ii] ]
			self.Ly = self.Ly + [ self.yhi[ii] - self.ylo[ii] ]
	
	def calc_displacement(self):
		for i in self.frames[1:]:
			self.df[i]['incUx']=self.df[i]["x"] -self.df[i-1]["x"]
			self.df[i]['cumUx']=self.df[i]["x"] -self.df[0]["x"]
			self.df[i]['incUy']=self.df[i]["y"] -self.df[i-1]["y"]
			self.df[i]['cumUy']=self.df[i]["y"] -self.df[0]["y"]

	def reg_mesh(self, mesh_size=2,trim_factor=0):
		self.mesh =[]
		for frame in self.frames[:-1]: 

			X = np.linspace(self.xlo[frame]+(trim_factor)*self.Lx[frame], 
				self.xhi[frame]-(trim_factor)*self.Ly[frame],
			 num=int(np.ceil(self.Lx[frame]*(1-2*trim_factor) /mesh_size)))

			Y = np.linspace(self.ylo[frame]+(trim_factor)*self.Ly[frame], 
				self.yhi[frame]-(trim_factor)*self.Ly[frame],
			 num=int(np.ceil(self.Ly[frame]*(1-2*trim_factor) /mesh_size)))

			self.mesh = self.mesh + [np.meshgrid(X,Y)] 

	def calc_strain_and_curl(self):

		self.calc_displacement()
		self.reg_mesh()


		incexx=[]
		inceyy=[]
		incexy=[]
		inceyx=[]
		self.inc_reg_Exx=[]
		self.inc_reg_Eyy=[]
		self.inc_reg_Exy=[]
		self.inc_reg_Ediv=[]
		self.inc_reg_curl=[]

		cumexx=[]
		cumeyy=[]
		cumexy=[]
		cumeyx=[]
		self.cum_reg_Exx=[]
		self.cum_reg_Eyy=[]
		self.cum_reg_Exy=[]
		self.cum_reg_Ediv=[]
		self.cum_reg_curl=[]
	    

		self.inc_reg_U=[]
		self.cum_reg_U=[]
		self.inc_par_Exx=[]
		self.inc_par_Eyy=[]
		self.inc_par_Exy=[]
		self.inc_par_Ediv=[]
		self.inc_par_curl=[]

		self.cum_par_Exx=[]
		self.cum_par_Eyy=[]
		self.cum_par_Exy=[]
		self.cum_par_Ediv=[]
		self.cum_par_curl=[]
	    

		for i in self.frames[1:]: 
			self.inc_reg_U = self.inc_reg_U + [[griddata((self.df[i]['x'].values, self.df[i]['y'].values),
	            self.df[i]['incUx'].values,
	            (self.mesh[i-1][0],self.mesh[i-1][1]), method='nearest'),griddata((self.df[i]['x'].values,
	            	self.df[i]['y'].values),self.df[i]['incUy'].values, 
	            (self.mesh[i-1][0],self.mesh[i-1][1]), method='nearest')]]
	                                                
			self.cum_reg_U = self.cum_reg_U + [[griddata((self.df[i]['x'].values, self.df[i]['y'].values),
	            self.df[i]['cumUx'].values,
	            (self.mesh[i-1][0],self.mesh[i-1][1]), method='nearest'),griddata((self.df[i]['x'].values,
	            	self.df[i]['y'].values),self.df[i]['cumUy'].values, 
	            (self.mesh[i-1][0],self.mesh[i-1][1]), method='nearest')]]                                
	                                                
		for i in self.frames[:-1]: 
			
			dexx, dexy = np.gradient(self.inc_reg_U[i][0])
			incexx = incexx+[dexx]
			incexy = incexy+[dexy]

			deyy, deyx = np.gradient(self.inc_reg_U[i][1])
			inceyy = inceyy+[deyy]
			inceyx = inceyx+[deyx]

			self.inc_reg_Exx = self.inc_reg_Exx + [incexx[i] + 0.5*(incexx[i]**2 + inceyx[i]**2)]                  
		        
			self.inc_reg_Eyy = self.inc_reg_Eyy + [inceyy[i] + 0.5*(inceyy[i]**2 + incexy[i]**2)]
		        
			self.inc_reg_Exy = self.inc_reg_Exy + [0.5*(incexy[i] + inceyx[i]) + 0.5*(incexx[i]*incexy[i] + inceyy[i]*inceyx[i])]
	        
			self.inc_reg_Ediv= self.inc_reg_Ediv + [np.sqrt((self.inc_reg_Exx[i]**2+self.inc_reg_Eyy[i]**2+self.inc_reg_Exy[i]**2+self.inc_reg_Exy[i]**2))]
	        
			self.inc_reg_curl = self.inc_reg_curl + [inceyx[i] - incexy[i]]



			dexx, dexy = np.gradient(self.cum_reg_U[i][0])
			cumexx = cumexx+[dexx]
			cumexy = cumexy+[dexy]

			deyy, deyx = np.gradient(self.cum_reg_U[i][1])
			cumeyy = cumeyy+[deyy]
			cumeyx = cumeyx+[deyx]

			self.cum_reg_Exx = self.cum_reg_Exx + [cumexx[i] + 0.5*(cumexx[i]**2 + cumeyx[i]**2)]                  
		        
			self.cum_reg_Eyy = self.cum_reg_Eyy + [cumeyy[i] + 0.5*(cumeyy[i]**2 + cumexy[i]**2)]
		        
			self.cum_reg_Exy = self.cum_reg_Exy + [0.5*(cumexy[i] + cumeyx[i]) + 0.5*(cumexx[i]*cumexy[i] + cumeyy[i]*cumeyx[i])]
	        
			self.cum_reg_Ediv= self.cum_reg_Ediv + [np.sqrt((self.cum_reg_Exx[i]**2+self.cum_reg_Eyy[i]**2+self.cum_reg_Exy[i]**2+self.cum_reg_Exy[i]**2))]
	        
			self.cum_reg_curl = self.cum_reg_curl + [cumeyx[i] - cumexy[i]]
	        	

		for i in self.frames[:-1]:
				
			temp = pd.DataFrame(np.column_stack([self.mesh[i][0].ravel(), self.mesh[i][1].ravel(), self.inc_reg_Exx[i].ravel() ]))
			self.inc_par_Exx =  self.inc_par_Exx + [ griddata((temp[0].values ,temp[1].values), temp[2].values, (self.df[i]['x'], self.df[i]['y']), method='nearest') ]
		        
			temp = pd.DataFrame(np.column_stack([self.mesh[i][0].ravel(), self.mesh[i][1].ravel(), self.inc_reg_Eyy[i].ravel() ]))
			self.inc_par_Eyy =  self.inc_par_Eyy + [ griddata((temp[0].values ,temp[1].values), temp[2].values, (self.df[i]['x'], self.df[i]['y']), method='nearest') ]
		        
			temp = pd.DataFrame(np.column_stack([self.mesh[i][0].ravel(), self.mesh[i][1].ravel(), self.inc_reg_Exy[i].ravel() ]))
			self.inc_par_Exy =  self.inc_par_Exy + [ griddata((temp[0].values ,temp[1].values), temp[2].values, (self.df[i]['x'], self.df[i]['y']), method='nearest') ]
		        
			temp = pd.DataFrame(np.column_stack([self.mesh[i][0].ravel(), self.mesh[i][1].ravel(), self.inc_reg_Ediv[i].ravel() ]))
			self.inc_par_Ediv =  self.inc_par_Ediv + [ griddata((temp[0].values ,temp[1].values), temp[2].values, (self.df[i]['x'], self.df[i]['y']), method='nearest') ]
		        
			temp = pd.DataFrame(np.column_stack([self.mesh[i][0].ravel(), self.mesh[i][1].ravel(), self.inc_reg_curl[i].ravel() ]))
			self.inc_par_curl =  self.inc_par_curl + [ griddata((temp[0].values ,temp[1].values), temp[2].values, (self.df[i]['x'], self.df[i]['y']), method='nearest') ]
		        

			temp = pd.DataFrame(np.column_stack([self.mesh[i][0].ravel(), self.mesh[i][1].ravel(), self.cum_reg_Exx[i].ravel() ]))
			self.cum_par_Exx =  self.cum_par_Exx + [ griddata((temp[0].values ,temp[1].values), temp[2].values, (self.df[i]['x'], self.df[i]['y']), method='nearest') ]
		        
			temp = pd.DataFrame(np.column_stack([self.mesh[i][0].ravel(), self.mesh[i][1].ravel(), self.cum_reg_Eyy[i].ravel() ]))
			self.cum_par_Eyy =  self.cum_par_Eyy + [ griddata((temp[0].values ,temp[1].values), temp[2].values, (self.df[i]['x'], self.df[i]['y']), method='nearest') ]
	        
			temp = pd.DataFrame(np.column_stack([self.mesh[i][0].ravel(), self.mesh[i][1].ravel(), self.cum_reg_Exy[i].ravel() ]))
			self.cum_par_Exy =  self.cum_par_Exy + [ griddata((temp[0].values ,temp[1].values), temp[2].values, (self.df[i]['x'], self.df[i]['y']), method='nearest') ]
		        
			temp = pd.DataFrame(np.column_stack([self.mesh[i][0].ravel(), self.mesh[i][1].ravel(), self.cum_reg_Ediv[i].ravel() ]))
			self.cum_par_Ediv =  self.cum_par_Ediv + [ griddata((temp[0].values ,temp[1].values), temp[2].values, (self.df[i]['x'], self.df[i]['y']), method='nearest') ]
		        
			temp = pd.DataFrame(np.column_stack([self.mesh[i][0].ravel(), self.mesh[i][1].ravel(), self.cum_reg_curl[i].ravel() ]))
			self.cum_par_curl =  self.cum_par_curl + [ griddata((temp[0].values ,temp[1].values), temp[2].values, (self.df[i]['x'], self.df[i]['y']), method='nearest') ]

class indentationYF:
	def __init__(self, xy_file, f_file, Ycolumns=['X', 'Y'], Fcolumns=['T','Fx', 'Fy']):
		
		self.xy = pd.read_csv(xy_file,skiprows=[0] ,sep='\s', names=Ycolumns,engine='python')

		self.f = pd.read_csv(f_file,skiprows=[0,1,2] ,sep='\s', names=Fcolumns,engine='python')

		self.xy['D'] = self.xy.loc[0,'Y'] - self.xy['Y'] 

def color_by_type(x, c1=(0.7,0.7,0.7), c2=(0.4,0.4,0.4)):
	if x ==4:
		return c1
	elif x ==3:
		return c1
	elif x==2:
		return c2
	else:
		return c2


def marker_sizer(xyz, r1=0.561231, r2tor1=1.5,i=0):
	fig = plt.gcf()
	size = fig.get_size_inches()
	Rho=3.14159265*round(xyz.N/2)*(1+r2tor1**2)*(r1**2)/(xyz.Lx[i]*xyz.Ly[i])
	markersize = Rho*size[0]*size[1]*(72**2)/(round(xyz.N/2)*(1+1.5**2))
	return markersize 		


def fig_sizer(xyz, scaler=100, i=0):

	w=xyz.Lx[i]/scaler
	l=xyz.Ly[i]/scaler

	return (w,l)


def show_mesh(xyz, reg=True, org=True, i=1, dpi=40, r1=0.561231, trim_factor=0.01, save=True):

	
	fig=plt.figure(i, figsize=fig_sizer(xyz,i=i), dpi=dpi)
	ax = fig.add_subplot(111)

	if org ==True:
		ax.scatter(xyz.df[i]['x'], xyz.df[i]['y'],s = marker_sizer(xyz, i=i, r1=r1)*xyz.df[i]['radius']/r1,
          	 c=xyz.df[i]['type'].apply(color_by_type).values)

	if reg ==True:
		ax.scatter(xyz.mesh[i][0], xyz.mesh[i][1], c='r')

	plt.axis('off')
	if save==True:
		plt.savefig('test.png', dpi=dpi)
	plt.show()


def plot_field(xyz, field, labels=['xlable', 'ylabel'], font=120, cmap='seismic',cbaro='Horizontal', org=True, axes=False, colorbar=False, 
	i=1, a1=0.5,a2=0.8, dpi=40, bmesh=True,  ws=False, r1=0.561231, save=True, save_file_name=None, save_directory=None,
	 legend=False, cmin=None, cmax=None):

	fig=plt.figure(figsize=fig_sizer(xyz, i=i), dpi=dpi)

	ax = fig.add_subplot(111)

	
	if ws == False:
            fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

	s = marker_sizer(xyz, 1)*xyz.df[1]['radius']/0.561231*0.38
	c = xyz.df[i]['type'].apply(color_by_type).values

	
	if bmesh==True:
            ax=plt.scatter(xyz.df[i].x, xyz.df[i].y, s= s, c = c, alpha=a2 )

	if (cmin == None) & (cmax == None):  
		cmin=np.mean(field[i])-3*np.std(field[i])
		cmax=np.mean(field[i])+3*np.std(field[i])   

	cax= plt.scatter(xyz.df[i]['x'], xyz.df[i]['y'], s = s,  c = field[i], vmin=cmin,
                         vmax=cmax, cmap=cmap, alpha =a1)

	if axes == False:
            plt.axis('off')
	if (colorbar == True):
		cbar = fig.colorbar(cax, orientation=cbaro,fraction=0.046, pad=0.04)
		cbar.ax.tick_params(labelsize=font/2)
		# plt.savefig(save_directory+'colorbar.png',dpi=dpi)

	plt.ylabel(labels[0], fontsize=font)
	plt.xlabel(labels[1], fontsize=font)
	plt.xticks(fontsize =font)
	plt.yticks(fontsize =font)

	if legend==True:
		plt.legend(loc=0, fontsize=font)
	
	
	if save == True:
		plt.savefig(save_directory+'\\'+save_file_name,dpi=dpi)

plt.show()


def show_field(xyz, field, labels=['xlable', 'ylabel'], org=True, i=1, dpi=40, r1=0.561231, save=True, save_file_name=None, save_directory=None):
	fig=plt.figure(figsize=fig_sizer(xyz, i=i), dpi=dpi)
	ax = fig.add_subplot(111)
	ax.imshow(field)
	plt.axis('off')

	plt.ylabel(labels[0], fontsize=120)
	plt.xlabel(labels[1], fontsize=120)
	plt.xticks(fontsize =120)
	plt.yticks(fontsize =120)
	plt.legend(loc=0, fontsize=100)
	if save==True:
		plt.savefig('teststrain.png', dpi=dpi)

	plt.show()