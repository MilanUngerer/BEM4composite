#Librerias
import numpy as np
import bempp.api

#########################################
#Parametros de trabajo (Input)
omega = 10. * 2.*np.pi * 1.e9
e0 = 8.854*1e-12 * 1e-18   #micrometros
mu0 = 4.*np.pi*1e-7 * 1e6  #micrometros
mu = (-18.-40.j) * mu0
cc = 6.5e4 * 1.e-18  #micrometros
antena = np.array([[1.e5, 5.e5, 10.e5],[0., 0., 0.],[0., 0., 0.]])
#antena = np.array([[1.e5],[0.],[0.]])
Amp = 1.e6  #micrometros
########################################


#calculo variables a utilizar
k0 = omega*np.sqrt(e0*mu0)
k = np.sqrt(1j*omega*mu*cc)
lam0 = 2*np.pi/k0
lam = 2*np.pi/k
alfa = mu/mu0

########################################

#Impresion de indices
print "Numero de onda exterior:", k0
print "Numero de onda interior conductor:", k
print "Indice de transmision conductor:", alfa
print "Longitud de onda interior:", lam, "micras"
print "Longitud de onda exterior:", lam0, "micras"

########################################

#Importando mallas
grid = bempp.api.import_grid('/home/milanungerer/meshes/45m_45mm/22624_45_45.msh')

########################################

#Funciones de dirichlet y neumann
def dirichlet_fun(x, n, domain_index, result):
	result[0] = Amp * np.exp(1j*k0*x[0])
def neumann_fun(x, n, domain_index, result):
	result[0] = Amp * 1j*k0*n[0]*np.exp(1j*k0*x[0])

########################################

#Operadores multitrazo
Ai = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid, k)
Ae = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid, k0)
ident = bempp.api.operators.boundary.sparse.multitrace_identity(grid)

#Transmision en Multitrazo
Ae[0,0] = Ae[0,0] * alfa
Ae[0,1] = Ae[0,1] * alfa
Ai[0,1] = Ai[0,1] * alfa
Ai[1,1] = Ai[1,1] * alfa

#Espacios
dirichlet_space = Ae[0,0].domain
neumann_space = Ae[0,1].domain

#LHS
blocked = 0.5 * ident * (alfa + 1) - Ai + Ae

#Condiciones de borde
dirichlet_grid_fun = bempp.api.GridFunction(dirichlet_space, fun=dirichlet_fun)
neumann_grid_fun = bempp.api.GridFunction(neumann_space, fun=neumann_fun)

#Discretizacion lado izquierdo
blocked_discretizado = blocked.strong_form()

#RHS
rhs = np.concatenate([alfa * dirichlet_grid_fun.coefficients, neumann_grid_fun.coefficients])

#Sistema de ecuaciones
import inspect
from scipy.sparse.linalg import gmres
array_it = np.array([])
array_frame = np.array([])
it_count = 0
def iteration_counter(x):
	global array_it
	global array_frame
	global it_count
	it_count += 1
	frame = inspect.currentframe().f_back
	array_it = np.append(array_it, it_count)
	array_frame = np.append(array_frame, frame.f_locals["resid"])
	if it_count % 10 == 0:
		print it_count, frame.f_locals["resid"]

print("Shape of matrix: {0}".format(blocked_discretizado.shape))
x,info = gmres(blocked_discretizado, rhs, tol=1e-5, callback = iteration_counter, maxiter = 5000, restart = 5000)
print("El sistema fue resuelto en {0} iteraciones".format(it_count))
np.savetxt("Solucion.out", x, delimiter=",")

#Campo exterior
exterior_field_dirichlet = bempp.api.GridFunction(dirichlet_space, coefficients=x[:dirichlet_space.global_dof_count])
exterior_field_neumann = bempp.api.GridFunction(neumann_space,coefficients=x[dirichlet_space.global_dof_count:])

#Campo interior
interior_field_dirichlet = exterior_field_dirichlet
interior_field_neumann = exterior_field_neumann * alfa

#Calculo campo en antena
slp_pot_ext = bempp.api.operators.potential.helmholtz.single_layer(dirichlet_space, antena, k0)
dlp_pot_ext = bempp.api.operators.potential.helmholtz.double_layer(dirichlet_space, antena, k0)

Campo_en_antena_dis = (dlp_pot_ext * exterior_field_dirichlet - slp_pot_ext * exterior_field_neumann).ravel() #+ Amp * np.exp(1j*k0*antena[0])
print "Valor del campo dispersado en receptor:", Campo_en_antena_dis
print "Valor del campo total en receptor:", Campo_en_antena_dis + Amp * np.exp(1j*k0*antena[0])

########### Graficos ######################

############################################################
#Grafico de convergencia
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from matplotlib import rcParams

#Grafico de convergencia
rcParams["font.family"] = "serif"
rcParams["font.size"] = 20
pyplot.figure(figsize = (15,10))
pyplot.title("Convergence")
pyplot.semilogy(array_it, array_frame, lw=2)
pyplot.xlabel("iteration")
pyplot.ylabel("residual")
pyplot.grid()
pyplot.savefig("Convergence.pdf")


#############################################################
import sys
sys.exit()
#############################################################


#Imagen con colores
Nx = 500
Ny = 500
xmin, xmax, ymin, ymax = [-1000, 1000, -1000, 1000]
plot_grid = np.mgrid[xmin:xmax:Nx * 1j, ymin:ymax:Ny * 1j]
points = np.vstack((plot_grid[0].ravel(), plot_grid[1].ravel(), np.zeros(plot_grid[0].size)))
u_evaluated = np.zeros(points.shape[1], dtype=np.complex128)

x, y = points[:2]
idx_ext = np.sqrt(x**2 + y**2) > 33.
idx_int = np.sqrt(x**2 + y**2) <= 33.

points_exterior = points[:, idx_ext]
points_interior = points[:, idx_int]

slp_pot_int = bempp.api.operators.potential.helmholtz.single_layer(dirichlet_space, points_interior, k)
slp_pot_ext = bempp.api.operators.potential.helmholtz.single_layer(dirichlet_space, points_exterior, k0)
dlp_pot_int = bempp.api.operators.potential.helmholtz.double_layer(dirichlet_space, points_interior, k)
dlp_pot_ext = bempp.api.operators.potential.helmholtz.double_layer(dirichlet_space, points_exterior, k0)

total_field_int = (slp_pot_int * interior_field_neumann - dlp_pot_int * interior_field_dirichlet).ravel()
total_field_ext = (dlp_pot_ext * exterior_field_dirichlet - slp_pot_ext * exterior_field_neumann).ravel() + Amp*np.exp(1j * k0 * points_exterior[0])

total_field = np.zeros(points.shape[1], dtype='complex128')
total_field[idx_ext] = total_field_ext
total_field[idx_int] = total_field_int
total_field = total_field.reshape([Nx, Ny])


##############################################################
import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab as plt

fig = plt.figure(figsize=(10, 8))
plt.imshow(np.real(total_field.T), extent=[-1000, 1000, -1000, 1000])
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.savefig('45m_perm1')


