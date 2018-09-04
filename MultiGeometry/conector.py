################################################################
###### ACA SE LEE EL PROGRAMA QUE CONTIENE LA INFORMACION ######
################################################################
import numpy as np

info = []
for line in open("info.txt"):
    li=line.strip()
    if not li.startswith("#"): #Comentarios empiezan con '#'
        info.append(line.split())

info = filter(None, info)
#print info
N_hilos = len(info) - 2
print '\nEl composito contiene', N_hilos, 'microhilos\n'


################################################################
###### ACA SE ESCRIBE EL PROGRAMA PARA EJECUTAR EL CALCULO #####
################################################################
exe = open('ejecutor.py','w')

###### IMPORTANDO LIBRERIAS ####################################
exe.write('#Preambulo\n')
exe.write('import numpy as np\n')
exe.write('import bempp.api\n')

###### PARAMETROS DE ENTRADA ###################################
#exe.write('bempp.api.global_parameters.quadrature.double_singular = 7\n') #Orden de cuadratura
exe.write('omega = 2.*np.pi*' + info[0][1] + '\n') #frecuencia angular
exe.write('e0 = 8.854*1e-12*1e-18\n') #permitividad del vacio
exe.write('mu0 = 4.*np.pi*1e-7*1e6\n') #permeabilidad del vacio
exe.write('mue = ' + info[0][2] + '*mu0\n') #permeabilidad de la matriz
exe.write('ee = ' + info[0][3] + '*e0\n') #permitividad de la matriz
exe.write('mui = ' + info[0][4] + '*mu0\n') #permeabilidad del conductor
exe.write('ei = ' + info[0][5] + '*e0\n') #permitividad del conductor
exe.write('k = omega*np.sqrt(e0*mu0)\n') #numero de onda exterior
exe.write('lam = 2*np.pi/k\n') #longitud de onda al exterior
exe.write('nm = np.sqrt((ee*mue)/(e0*mu0))\n') #indice de refraccion matriz
exe.write('nc = np.sqrt((ei*mui)/(e0*mu0))\n') #indice de refraccion conductor
exe.write('alfa_m = mue/mu0\n') #indice de transmision matriz
exe.write('alfa_c = mui/mue\n') #indice de transmision conductor
exe.write('antena = np.array('+info[0][6]+')\n') #Punto antena receptor

###### ESCRIBIR VALORES DE INTERES #############################
exe.write('print "Numero de onda exterior:", k\n')
exe.write('print "Indice de refraccion matriz:", nm\n')
exe.write('print "Indice de refraccion conductor:", nc\n')
exe.write('print "Numero de onda interior matriz:", nm*k\n')
exe.write('print "Numero de onda interior conductor:", nm*nc*k\n')
exe.write('print "Indice de transmision matriz:", alfa_m\n')
exe.write('print "Indice de transmision conductor:", alfa_c\n')
exe.write('print "Longitud de onda:", lam, "micras"\n')

###### IMPORTANDO MALLAS #######################################
exe.write('\n#Importando mallas\n')

exe.write('matriz = bempp.api.import_grid('+info[1][0]+')\n') #Malla de la matriz
for i in range(N_hilos): #Mallas de hilos
    exe.write('grid_'+str(i)+' = bempp.api.import_grid('+info[i+2][0]+')\n') 

###### FUNCIONES DIRICHLET Y NEUMANN ###########################
exe.write('\n#Funciones de dirichlet y neumann\n')

exe.write('def dirichlet_fun(x, n, domain_index, result):\n')
exe.write('\tresult[0] = '+info[0][0]+'*np.exp(1j*k*x[0])\n')

exe.write('def neumann_fun(x, n, domain_index, result):\n')
exe.write('\tresult[0] = '+info[0][0]+'*1j*k*n[0]*np.exp(1j*k*x[0])\n')

###### OPERADORES EN EL BORDE ##################################
exe.write('\n#Operadores multitrazo\n') 
exe.write('Ai_m = bempp.api.operators.boundary.helmholtz.multitrace_operator(matriz, nm*k)\n') #Operador multitrazo interior matriz
exe.write('Ae_m = bempp.api.operators.boundary.helmholtz.multitrace_operator(matriz, k)\n') #Operador multitrazo exterior matriz
for i in range(N_hilos): #Operadores multitrazo hilos
    exe.write('Ai_'+str(i)+' = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid_' + str(i) + ',nm*nc*k)\n') #Interior hilos
    exe.write('Ae_'+str(i)+' = bempp.api.operators.boundary.helmholtz.multitrace_operator(grid_' + str(i) + ',nm*k)\n') #Exterior hilos

exe.write('\n#Transmision en Multitrazo\n')

exe.write('Ae_m[0,1] = Ae_m[0,1]*(1./alfa_m)\n') #Transmision en matriz exterior
exe.write('Ae_m[1,1] = Ae_m[1,1]*(1./alfa_m)\n') #Transmision en matriz exterior
for i in range(N_hilos): #Transmision en Multitrazo en hilos
    exe.write('Ai_'+str(i)+'[0,1] = Ai_'+str(i)+'[0,1]*alfa_c\n') #Transmision interior hilos
    exe.write('Ai_'+str(i)+'[1,1] = Ai_'+str(i)+'[1,1]*alfa_c\n') #Transmision interior hilos

exe.write('\n#Acople interior y exterior\n')
exe.write('op_m = (Ai_m + Ae_m)\n') #Interior + exterior matriz
for i in range(N_hilos): #Interior + exterior hilos
    exe.write('op_'+str(i)+' = (Ai_'+str(i)+' + Ae_'+str(i)+')\n')

exe.write('\n#Espacios\n')
exe.write('dirichlet_space_m = Ai_m[0,0].domain\n') #Espacio de dirichlet en matriz
exe.write('neumann_space_m = Ai_m[0,1].domain\n') #Espacio de neumann en matriz
for i in range(N_hilos): #Espacios en hilos
    exe.write('dirichlet_space_'+str(i)+' = Ai_'+str(i)+'[0,0].domain\n') #Espacio de dirichlet en hilos
    exe.write('neumann_space_'+str(i)+' = Ai_'+str(i)+'[0,1].domain\n') #Espacio de neumann en hilos

#operadores identidad
exe.write('\n#Operadores identidad\n')
exe.write('ident_m = bempp.api.operators.boundary.sparse.identity(neumann_space_m, neumann_space_m, neumann_space_m)\n') #Operador identidad matriz
for i in range(N_hilos): #Operadores identidad en hilos
    exe.write('ident_'+str(i)+' = bempp.api.operators.boundary.sparse.identity(neumann_space_' + str(i) + ', neumann_space_' + str(i) + ', neumann_space_' + str(i) + ')\n') #Identidad

#operadores diagonales
exe.write('\n#Operadores diagonales\n')
exe.write('op_m[1,1] = op_m[1,1] + 0.5 * ident_m * ((alfa_m -1)/alfa_m)\n')
for i in range(N_hilos):
    exe.write('op_'+str(i)+'[1,1] = op_' + str(i) + '[1,1] + 0.5 * ident_' + str(i) + '* (alfa_c - 1)\n')

#Operadores entre mallas
exe.write('\n#Operadores entre mallas\n')
for i in range(N_hilos): #Operadores entre mallas 
    exe.write('SLP_m_'+str(i)+' = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_m, dirichlet_space_'+str(i)+', dirichlet_space_'+str(i)+', nm*k)\n') #Operadores matriz-hilos single layer
    exe.write('SLP_'+str(i)+'_m = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_'+str(i)+', dirichlet_space_m, dirichlet_space_m, nm*k)\n') #Operadores matriz-hilos single layer

    exe.write('DLP_m_'+str(i)+' = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_m, dirichlet_space_'+str(i)+', dirichlet_space_'+str(i)+', nm*k)\n') #Operadores matriz-hilos double layer
    exe.write('DLP_'+str(i)+'_m = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_'+str(i)+', dirichlet_space_m, dirichlet_space_m, nm*k)\n') #Operadores matriz-hilos double layer

    exe.write('ADLP_m_'+str(i)+' = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_m, neumann_space_'+str(i)+', neumann_space_'+str(i)+', nm*k)\n') #Operadores matriz-hilos adjoint double layer
    exe.write('ADLP_'+str(i)+'_m = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_'+str(i)+', neumann_space_m, neumann_space_m, nm*k)\n') #Operadores matriz-hilos adjoint double layer

    exe.write('HYP_m_'+str(i)+' = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_m, neumann_space_'+str(i)+', neumann_space_'+str(i)+', nm*k)\n') #Operadores matriz-hilos hypersingular
    exe.write('HYP_'+str(i)+'_m = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_'+str(i)+', neumann_space_m, neumann_space_m, nm*k)\n') #Operadores matriz-hilos hypersingular

    for j in range(N_hilos): #Interaccion entre hilos
        if i!=j:
            exe.write('SLP_'+str(i)+'_'+str(j)+' = bempp.api.operators.boundary.helmholtz.single_layer(neumann_space_'+str(i)+', dirichlet_space_'+str(j)+', dirichlet_space_'+str(j)+', nm*k)\n') #Single-layer interaccion entre hilos
            exe.write('DLP_'+str(i)+'_'+str(j)+' = bempp.api.operators.boundary.helmholtz.double_layer(dirichlet_space_'+str(i)+', dirichlet_space_'+str(j)+', dirichlet_space_'+str(j)+', nm*k)\n') #Double-layer interaccion entre hilos
            exe.write('ADLP_'+str(i)+'_'+str(j)+' = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(neumann_space_'+str(i)+', neumann_space_'+str(j)+', neumann_space_'+str(j)+', nm*k)\n') #Adjoint interaccion entre hilos
            exe.write('HYP_'+str(i)+'_'+str(j)+' = bempp.api.operators.boundary.helmholtz.hypersingular(dirichlet_space_'+str(i)+', neumann_space_'+str(j)+', neumann_space_'+str(j)+', nm*k)\n') #Hypersingular interaccion entre hilos


###### ENSAMBLANDO MATRIZ DE OPERADORES ########################
exe.write('\n#Matriz de operadores\n')
exe.write('blocked = bempp.api.BlockedOperator('+str(2*(N_hilos+1))+','+str(2*(N_hilos+1))+')\n') #Tamano bloque de operadores

exe.write('\n#Diagonal\n')
exe.write('blocked[0,0] = op_m[0,0]\n') #Diagonal matriz
exe.write('blocked[0,1] = op_m[0,1]\n') #Diagonal matriz
exe.write('blocked[1,0] = op_m[1,0]\n') #Diagonal matriz
exe.write('blocked[1,1] = op_m[1,1]\n') #Diagonal matriz

c=0
for i in range(2, 2*(N_hilos+1)-1, 2): #Diagonal hilos
    exe.write('blocked['+str(i)+','+str(i)+'] = op_'+str(c)+'[0,0]\n')
    exe.write('blocked['+str(i)+','+str(i+1)+'] = op_'+str(c)+'[0,1]\n')
    exe.write('blocked['+str(i+1)+','+str(i)+'] = op_'+str(c)+'[1,0]\n') 
    exe.write('blocked['+str(i+1)+','+str(i+1)+'] = op_'+str(c)+'[1,1]\n')
    c+=1
        
exe.write('\n#Contribucion hilos-matriz\n')
c=0
for i in range(2, 2*N_hilos+1, 2): #Contribucion hilos en matriz
    exe.write('blocked[0,'+str(i)+'] = DLP_'+str(c)+'_m\n') #Double-layer hilos-matriz 
    exe.write('blocked[0,'+str(i+1)+'] = -SLP_'+str(c)+'_m\n') #Single-layer hilos-matriz
    exe.write('blocked[1,'+str(i)+'] = -HYP_'+str(c)+'_m\n') #Hypersingular hilos-matriz
    exe.write('blocked[1,'+str(i+1)+'] = -ADLP_'+str(c)+'_m\n') #Adjoint hilos-matriz
    c+=1
        
c1=0
for i in range(2, 2*(N_hilos+1)-1, 2): #Contribucion hilos-hilos
    c2=0
    for j in range(2, 2*(N_hilos+1)-1, 2):
        if i<j:
            exe.write('\n#Contribucion hilos-hilos\n')
            exe.write('blocked['+str(i)+','+str(j)+'] = DLP_'+str(c2+1)+'_'+str(c1)+'\n') #Double-layer hilo-hilo
            exe.write('blocked['+str(i)+','+str(j+1)+'] = -SLP_'+str(c2+1)+'_'+str(c1)+'\n') #Single-layer hilo-hilo
            exe.write('blocked['+str(i+1)+','+str(j)+'] = -HYP_'+str(c2+1)+'_'+str(c1)+'\n') #Hypersingular hilo-hilo
            exe.write('blocked['+str(i+1)+','+str(j+1)+'] = -ADLP_'+str(c2+1)+'_'+str(c1)+'\n') #Adjoint hilo-hilo
            c2+=1
        elif i>j:
            exe.write('\n#Contribucion hilos-hilos\n')
            exe.write('blocked['+str(i)+','+str(j)+'] = DLP_'+str(c2)+'_'+str(c1)+'\n') #Double-layer hilo-hilo
            exe.write('blocked['+str(i)+','+str(j+1)+'] = -SLP_'+str(c2)+'_'+str(c1)+'\n') #Single-layer hilo-hilo
            exe.write('blocked['+str(i+1)+','+str(j)+'] = -HYP_'+str(c2)+'_'+str(c1)+'\n') #Hypersingular hilo-hilo
            exe.write('blocked['+str(i+1)+','+str(j+1)+'] = -ADLP_'+str(c2)+'_'+str(c1)+'\n') #Adjoint hilo-hilo
            c2+=1
    
    exe.write('\n#Contribucion matriz-hilos\n') #Contribucion matriz en hilos
    exe.write('blocked['+str(i)+',0] = -DLP_m_'+str(c1)+'\n') #Double-layer matriz-hilos
    exe.write('blocked['+str(i)+',1] = SLP_m_'+str(c1)+'\n') #Single-layer matriz-hilos
    exe.write('blocked['+str(i+1)+',0] = HYP_m_'+str(c1)+'\n') #Hypersingular matriz-hilos 
    exe.write('blocked['+str(i+1)+',1] = ADLP_m_'+str(c1)+'\n') #Adjoint matriz-hilos
    c1+=1

###### CONDICIONES DE BORDE ####################################
exe.write('\n#Condiciones de borde\n') #Condiciones en el borde de la matriz
exe.write('dirichlet_grid_fun_m = bempp.api.GridFunction(dirichlet_space_m, fun=dirichlet_fun)\n')
exe.write('neumann_grid_fun_m = bempp.api.GridFunction(neumann_space_m, fun=neumann_fun)\n')

###### DISCRETIZACION ##########################################
exe.write('\n#Discretizacion lado izquierdo\n') #Lado izquierdo
exe.write('blocked_discretizado = blocked.strong_form()\n')

exe.write('\n#Discretizacion lado derecho\n') #Lado derecho con onda incidente
exe.write('rhs = np.concatenate([')
exe.write('dirichlet_grid_fun_m.coefficients, neumann_grid_fun_m.coefficients,')
for i in range(N_hilos):
    exe.write('np.zeros(dirichlet_space_'+str(i)+'.global_dof_count), np.zeros(neumann_space_'+str(i)+'.global_dof_count)')
    if i!=N_hilos-1:
        exe.write(', ')
exe.write('])\n')


###### SISTEMA DE ECUACIONES ###################################
exe.write('\n#Sistema de ecuaciones\n')
exe.write('import inspect\n')
exe.write('from scipy.sparse.linalg import gmres\n')

exe.write('array_it = np.array([])\n')
exe.write('array_frame = np.array([])\n')
exe.write('it_count = 0\n') #numero de iteraciones
exe.write('def iteration_counter(x):\n') #Contador
exe.write('\tglobal array_it\n')
exe.write('\tglobal array_frame\n')
exe.write('\tglobal it_count\n')
exe.write('\tit_count += 1\n')
exe.write('\tframe = inspect.currentframe().f_back\n')
exe.write('\tarray_it = np.append(array_it, it_count)\n')
exe.write('\tarray_frame = np.append(array_frame, frame.f_locals["resid"])\n')
#exe.write('\tprint it_count, frame.f_locals["resid"]\n') #Descomentar para ver numero iteracion y convergencia en cada iteracion

exe.write('print("Shape of matrix: {0}".format(blocked_discretizado.shape))\n') #Tamano de la matriz
exe.write('x,info = gmres(blocked_discretizado, rhs, tol=1e-5, callback = iteration_counter, maxiter = 50000, restart = 5000)\n') #GMRES para resolver el sistema lineal
exe.write('print("El sistema fue resuelto en {0} iteraciones".format(it_count))\n') #Numero de iteraciones

exe.write('np.savetxt("Solucion.out", x, delimiter=",")\n') #Guardando solucion en archivo txt

###### SEPARACION DE LA SOLUCION ###############################
exe.write('\n#Campo interior\n') #Separar la solucion del sistema solo para la matriz
exe.write('interior_field_dirichlet_m = bempp.api.GridFunction(dirichlet_space_m, coefficients=x[:dirichlet_space_m.global_dof_count])\n')
exe.write('interior_field_neumann_m = bempp.api.GridFunction(neumann_space_m,coefficients=x[dirichlet_space_m.global_dof_count:dirichlet_space_m.global_dof_count + neumann_space_m.global_dof_count])\n')

###### CAMPO EXTERIOR ##########################################
exe.write('\n#Campo exterior\n') #Aplicar condiciones de transmision para obtener campo exterior 
exe.write('exterior_field_dirichlet_m = interior_field_dirichlet_m\n')
exe.write('exterior_field_neumann_m = interior_field_neumann_m*(1./alfa_m)\n')

###### CALCULO DEL CAMPO EN ANTENA #############################
exe.write('\n#Calculo campo en antena\n')

###### CAMPO EXTERIOR A LA MATRIZ ##############################
exe.write('slp_pot_ext_m = bempp.api.operators.potential.helmholtz.single_layer(dirichlet_space_m, antena, k)\n') #Single-layer exterior
exe.write('dlp_pot_ext_m = bempp.api.operators.potential.helmholtz.double_layer(dirichlet_space_m, antena, k)\n') #Double-layer exterior

exe.write('Campo_en_antena = (dlp_pot_ext_m * exterior_field_dirichlet_m - slp_pot_ext_m * exterior_field_neumann_m).ravel() + np.exp(1j*k*antena[0])\n') #Calculo del campo en la antena

exe.write('print "Valor del campo en receptor:", Campo_en_antena\n') #Imprimir resultados

###### GRAFICO DE CONVERGENCIA #################################
exe.write('\nimport matplotlib\n')
exe.write('matplotlib.use("Agg")\n')
exe.write('from matplotlib import pyplot\n')
exe.write('from matplotlib import rcParams\n')
exe.write('rcParams["font.family"] = "serif"\n')
exe.write('rcParams["font.size"] = 20\n')
exe.write('pyplot.figure(figsize = (15,10))\n')
exe.write('pyplot.title("Convergence")\n')
exe.write('pyplot.plot(array_it, array_frame, lw=2)\n')
exe.write('pyplot.xlabel("iteration")\n')
exe.write('pyplot.ylabel("residual")\n')
exe.write('pyplot.grid()\n')
exe.write('pyplot.savefig("Convergence.pdf")\n')

###### CERRAMOS ARCHIVO .PY CREADO #############################
exe.close

################################################################
