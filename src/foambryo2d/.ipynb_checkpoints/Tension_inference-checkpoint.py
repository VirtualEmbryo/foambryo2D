import numpy as np 
from scipy import linalg
import matplotlib.pyplot as plt


def change_idx_labels(labels,mapping):
    #mapping : {key_init:key_end} has to be a bijection
    new_faces = labels.copy()
    for key in mapping : 
        new_faces[labels==key]=mapping[key]
    return(new_faces)

def change_idx_cells(faces,mapping):
    #mapping : {key_init:key_end} has to be a bijection
    new_faces = faces.copy()
    for key in mapping : 
        new_faces[faces[:,3]==key][x:,3]=mapping[key]
        new_faces[faces[:,4]==key][:,4]=mapping[key]
    return(new_faces)


def infer_tension(Mesh,mean_tension=1,check_topo_changes=False,mode='YD'):

    if mode == 'YD': 
        return(infer_tension_symmetrical_yd(Mesh,mean_tension,check_topo_changes))
    elif mode == 'Projection_YD': 
        return(infer_tension_projection_yd(Mesh,mean_tension,check_topo_changes))
    elif mode == 'cotan': 
        return(infer_tension_cotan(Mesh,mean_tension,check_topo_changes))
    elif mode == 'inv_cotan':
        return(infer_tension_inv_cotan(Mesh,mean_tension,check_topo_changes))
    elif mode =='Lamy': 
        return(infer_tension_lamy(Mesh,mean_tension,check_topo_changes))
    elif mode =='inv_Lamy': 
        return(infer_tension_inv_lamy(Mesh,mean_tension,check_topo_changes))
    elif mode =='Lamy_Log': 
        return(infer_tension_lamy_log(Mesh,mean_tension,check_topo_changes))
    else : 
        print("Unimplemented method")
        return(None)

    
def check_trijunctions(dict_rad,dict_lengths,dict_tensions):
    print("Check of topology changes...")
        
    T = np.array(sorted(list(dict_lengths.keys())))
    kt = np.amax(T)+1
    Keys_t = T[:,0]*kt + T[:,1]
    reverse_map_t = dict(zip(Keys_t, np.arange(len(Keys_t))))
    
    A = np.array(list(dict_rad.keys()))
    A = -np.sort(-A,axis=1)
    ka = np.amax(A)
    Keys_a = A[:,0]*ka**2 + A[:,1]*ka + A[:,2]
    Keys_a,index = np.unique(Keys_a,return_index=True)
    A=A[index]
    
    list_topo_change = []
    for i,key in enumerate(Keys_a):
        c = key//ka**2
        b = (key-(c*(ka**2)))//ka
        a = key%ka
        
        idx_tab = reverse_map_t[min(a,b)*kt+max(a,b)]
        idx_tbc = reverse_map_t[min(b,c)*kt+max(b,c)]
        idx_tac = reverse_map_t[min(a,c)*kt+max(a,c)]
        
        tab = dict_tensions[(min(a,b),max(a,b))]
        tbc = dict_tensions[(min(b,c),max(b,c))]
        tac = dict_tensions[(min(a,c),max(a,c))]
        
        if not(tab+tbc-tac>0 and tab+tac-tbc>0 and tbc+tac-tab>0) : 
            list_topo_change.append([a,b,c])
            
    if len(list_topo_change)>0 : 
            print("Problematic junctions:",list_topo_change)
    else : 
        print("No topology changes detected !")
    return(list_topo_change)
        


def infer_forces(Mesh, mean_tension=1, P0 = 0, mode_tension = "YD", mode_pressure = "Laplace"): 
    _,dict_tensions,_ = infer_tension(Mesh, mean_tension=mean_tension, mode = mode_tension)
    _,dict_pressures,_ = infer_pressure(Mesh,dict_tensions,mode=mode_pressure, P0 = P0)
    return(dict_tensions, dict_pressures)

def change_idx_cells(faces,mapping):
    #mapping : {key_init:key_end} has to be a bijection
    new_faces = faces.copy()
    for key in mapping : 
        new_faces[faces[:,3]==key][x:,3]=mapping[key]
        new_faces[faces[:,4]==key][:,4]=mapping[key]
    return(new_faces)



def cot(x): 
    return(np.cos(x)/np.sin(x))
def compute_z(phi1,phi2,phi3): 
    z1 = np.sin(phi1)/(1-np.cos(phi1))
    z2 = np.sin(phi2)/(1-np.cos(phi2))
    z3 = np.sin(phi3)/(1-np.cos(phi3))
    return(z1,z2,z3)

def factors_z(z1,z2,z3): 
    f1 = z2*z3
    f2 = z1*z3
    f3 = z1*z2
    s1 = z1*(z2 + z3)/2
    s2 = z2*(z1 + z3)/2
    s3 = z3*(z1 + z2)/2
    return(f1,f2,f3,s1,s2,s3)






def build_matrix_tension_symmetrical_yd(dict_rad,dict_lengths,n_materials,mean_tension=1):

    #Index interfaces :
    T = np.array(sorted(list(dict_lengths.keys())))
    kt = np.amax(T)+1
    Keys_t = T[:,0]*kt + T[:,1]
    reverse_map_t = dict(zip(Keys_t, np.arange(len(Keys_t))))
    
    #Index angles : 
    A = np.array(list(dict_rad.keys()))
    A = -np.sort(-A,axis=1)
    ka = np.amax(A)
    Keys_a = A[:,0]*ka**2 + A[:,1]*ka + A[:,2]
    Keys_a,index = np.unique(Keys_a,return_index=True)
    A=A[index]
    
    #MX = B
    #x : structure : [0,nm[ : tensions
    nm = len(T)
    nc = n_materials-1
    nj = len(A)
    size = (3*nj + 1)*(nm)
    B = np.zeros(3*nj + 1)
    B[-1]=nm*mean_tension  #Sum of tensions = number of membrane <=> mean of tensions = 1
    M = np.zeros((3*nj + 1, nm))
    
    
    #CLASSICAL 
    #YOUNG-DUPRÉ
    for i,key in enumerate(Keys_a):
        c = key//ka**2
        b = (key-(c*(ka**2)))//ka
        a = key%ka

        O_b = dict_rad[(a,b,c)] #O_b = abc angle
        O_a = dict_rad[(b,a,c)] #O_a = cab angle
        O_c = dict_rad[(a,c,b)] #O_c = acb angle
        
        idx_tab = reverse_map_t[min(a,b)*kt+max(a,b)]
        idx_tbc = reverse_map_t[min(b,c)*kt+max(b,c)]
        idx_tac = reverse_map_t[min(a,c)*kt+max(a,c)]

        M[3*i,idx_tab]=1             
        M[3*i,idx_tbc]=np.cos(O_b)  
        M[3*i,idx_tac]=np.cos(O_a)   

        M[3*i+1,idx_tab]=np.cos(O_b)  
        M[3*i+1,idx_tbc]=1            
        M[3*i+1,idx_tac]=np.cos(O_c)  

        M[3*i+2,idx_tab]=np.cos(O_a)            
        M[3*i+2,idx_tbc]=np.cos(O_c) 
        M[3*i+2,idx_tac]=1               

    M[3*nj,:nm]=1 #Enforces mean of tensions = 1 (see B[-1])

    return(M,B)



def build_matrix_tension_inv_lamy(dict_rad,dict_lengths,n_materials,mean_tension=1):

    #Index interfaces :
    T = np.array(sorted(list(dict_lengths.keys())))
    kt = np.amax(T)+1
    Keys_t = T[:,0]*kt + T[:,1]
    reverse_map_t = dict(zip(Keys_t, np.arange(len(Keys_t))))
    
    #Index angles : 
    A = np.array(list(dict_rad.keys()))
    A = -np.sort(-A,axis=1)
    ka = np.amax(A)
    Keys_a = A[:,0]*ka**2 + A[:,1]*ka + A[:,2]
    Keys_a,index = np.unique(Keys_a,return_index=True)
    A=A[index]
    
    #MX = B
    #x : structure : [0,nm[ : tensions
    nm = len(T)
    nc = n_materials-1
    nj = len(A)
    size = (3*nj + 1)*(nm)
    B = np.zeros(3*nj + 1)
    B[-1]=nm*mean_tension  #Sum of tensions = number of membrane <=> mean of tensions = 1
    M = np.zeros((3*nj + 1, nm))
    
    
    #CLASSICAL 
    #YOUNG-DUPRÉ
    for i,key in enumerate(Keys_a):
        c = key//ka**2
        b = (key-(c*(ka**2)))//ka
        a = key%ka

        O_b = dict_rad[(a,b,c)] #O_b = abc angle
        O_a = dict_rad[(b,a,c)] #O_a = cab angle
        O_c = dict_rad[(a,c,b)]
        idx_tab = reverse_map_t[min(a,b)*kt+max(a,b)]
        idx_tbc = reverse_map_t[min(b,c)*kt+max(b,c)]
        idx_tac = reverse_map_t[min(a,c)*kt+max(a,c)]

        M[3*i,idx_tab]=np.sin(O_a)    
        M[3*i,idx_tbc]=-np.sin(O_c)   

        M[3*i+1,idx_tab]=np.sin(O_b)  
        M[3*i+1,idx_tac]=-np.sin(O_c) 

        M[3*i+2,idx_tbc]=np.sin(O_b)  
        M[3*i+2,idx_tac]=-np.sin(O_a) 

    M[3*nj,:nm]=1 #Enforces mean of tensions = 1 (see B[-1])

    return(M,B)





def build_matrix_tension_lamy(dict_rad,dict_lengths,n_materials,mean_tension=1):

    #Index interfaces :
    T = np.array(sorted(list(dict_lengths.keys())))
    kt = np.amax(T)+1
    Keys_t = T[:,0]*kt + T[:,1]
    reverse_map_t = dict(zip(Keys_t, np.arange(len(Keys_t))))
    
    #Index angles : 
    A = np.array(list(dict_rad.keys()))
    A = -np.sort(-A,axis=1)
    ka = np.amax(A)
    Keys_a = A[:,0]*ka**2 + A[:,1]*ka + A[:,2]
    Keys_a,index = np.unique(Keys_a,return_index=True)
    A=A[index]
    
    #MX = B
    #x : structure : [0,nm[ : tensions
    nm = len(T)
    nc = n_materials-1
    nj = len(A)
    size = (3*nj + 1)*(nm)
    B = np.zeros(3*nj + 1)
    B[-1]=nm*mean_tension  #Sum of tensions = number of membrane <=> mean of tensions = 1
    M = np.zeros((3*nj + 1, nm))
    
    
    #CLASSICAL 
    #YOUNG-DUPRÉ
    for i,key in enumerate(Keys_a):
        c = key//ka**2
        b = (key-(c*(ka**2)))//ka
        a = key%ka

        O_b = dict_rad[(a,b,c)] #O_b = abc angle
        O_a = dict_rad[(b,a,c)] #O_a = cab angle
        O_c = dict_rad[(a,c,b)]
        idx_tab = reverse_map_t[min(a,b)*kt+max(a,b)]
        idx_tbc = reverse_map_t[min(b,c)*kt+max(b,c)]
        idx_tac = reverse_map_t[min(a,c)*kt+max(a,c)]

        M[3*i,idx_tab]=1/np.sin(O_c)    
        M[3*i,idx_tbc]=-1/np.sin(O_a)   

        M[3*i+1,idx_tab]=1/np.sin(O_c)  
        M[3*i+1,idx_tac]=-1/np.sin(O_b) 

        M[3*i+2,idx_tbc]=1/np.sin(O_a)  
        M[3*i+2,idx_tac]=-1/np.sin(O_b) 

    M[3*nj,:nm]=1 #Enforces mean of tensions = 1 (see B[-1])

    return(M,B)


def build_matrix_tension_projection_yd(dict_rad,dict_lengths,n_materials,mean_tension=1):

    #Index interfaces :
    T = np.array(sorted(list(dict_lengths.keys())))
    kt = np.amax(T)+1
    Keys_t = T[:,0]*kt + T[:,1]
    reverse_map_t = dict(zip(Keys_t, np.arange(len(Keys_t))))
    
    #Index angles : 
    A = np.array(list(dict_rad.keys()))
    A = -np.sort(-A,axis=1)
    ka = np.amax(A)
    Keys_a = A[:,0]*ka**2 + A[:,1]*ka + A[:,2]
    Keys_a,index = np.unique(Keys_a,return_index=True)
    A=A[index]
    
    #MX = B
    #x : structure : [0,nm[ : tensions
    nm = len(T)
    nc = n_materials-1
    nj = len(A)
    size = (2*nj + 1)*(nm)
    B = np.zeros(2*nj + 1)
    B[-1]=nm*mean_tension  #Sum of tensions = number of membrane <=> mean of tensions = 1
    M = np.zeros((2*nj + 1, nm))
    
    
    #CLASSICAL 
    #YOUNG-DUPRÉ
    for i,key in enumerate(Keys_a):
        c = key//ka**2
        b = (key-(c*(ka**2)))//ka
        a = key%ka
        O_b = dict_rad[(a,b,c)] #O_b = abc angle
        O_a = dict_rad[(b,a,c)] #O_a = cab angle
        idx_tab = reverse_map_t[min(a,b)*kt+max(a,b)]
        idx_tbc = reverse_map_t[min(b,c)*kt+max(b,c)]
        idx_tac = reverse_map_t[min(a,c)*kt+max(a,c)]

        M[2*i,idx_tab]=1              
        M[2*i,idx_tbc]=np.cos(O_b)    
        M[2*i,idx_tac]=np.cos(O_a)    
        M[2*i+1,idx_tac]=-np.sin(O_a) 
        M[2*i+1,idx_tbc]=np.sin(O_b)  

    M[2*nj,:nm]=1 #Enforces mean of tensions = 1 (see B[-1])

    return(M,B)









def build_matrix_tension_lamy_log(dict_rad,dict_lengths,n_materials,mean_tension=1):

    #Index interfaces :
    T = np.array(sorted(list(dict_lengths.keys())))
    kt = np.amax(T)+1
    Keys_t = T[:,0]*kt + T[:,1]
    reverse_map_t = dict(zip(Keys_t, np.arange(len(Keys_t))))
    
    #Index angles : 
    A = np.array(list(dict_rad.keys()))
    A = -np.sort(-A,axis=1)
    ka = np.amax(A)
    Keys_a = A[:,0]*ka**2 + A[:,1]*ka + A[:,2]
    Keys_a,index = np.unique(Keys_a,return_index=True)
    A=A[index]
    
    #MX = B
    #x : structure : [0,nm[ : tensions
    nm = len(T)
    nc = n_materials-1
    nj = len(A)
    size = (3*nj + 1)*(nm)
    B = np.zeros(3*nj + 1)
    B[-1]=nm*np.log(mean_tension)  #Sum of tensions = number of membrane <=> mean of tensions = 1
    M = np.zeros((3*nj + 1, nm))
    
    
    #CLASSICAL 
    #YOUNG-DUPRÉ
    for i,key in enumerate(Keys_a):
        c = key//ka**2
        b = (key-(c*(ka**2)))//ka
        a = key%ka

        O_b = dict_rad[(a,b,c)] #O_b = abc angle
        O_a = dict_rad[(b,a,c)] #O_a = cab angle
        O_c = dict_rad[(a,c,b)]
        idx_tab = reverse_map_t[min(a,b)*kt+max(a,b)]
        idx_tbc = reverse_map_t[min(b,c)*kt+max(b,c)]
        idx_tac = reverse_map_t[min(a,c)*kt+max(a,c)]

        M[3*i,idx_tab]=1    
        M[3*i,idx_tbc]=-1   
        B[3*i] = (np.log(np.sin(O_c)) - np.log(np.sin(O_a)))

        M[3*i+1,idx_tab]=1    
        M[3*i+1,idx_tac]=-1   
        B[3*i+1] = (np.log(np.sin(O_c)) - np.log(np.sin(O_b)))

        M[3*i+2,idx_tbc]=1    
        M[3*i+2,idx_tac]=-1   
        B[3*i+2] = (np.log(np.sin(O_a)) - np.log(np.sin(O_b)))

    M[3*nj,:nm]=1 #Enforces mean of tensions = 1 (see B[-1])

    return(M,B)


def cot(x): 
    return(np.cos(x)/np.sin(x))
def tensions_equilibrium(phi1,phi2,phi3):
    t1 = cot(phi1/2)*(cot(phi2/2) + cot(phi3/2))
    t2 = cot(phi2/2)*(cot(phi1/2) + cot(phi3/2))
    t3 = cot(phi3/2)*(cot(phi1/2) + cot(phi2/2))
    return(t1,t2,t3)


def build_matrix_tension_cotan(dict_rad,dict_lengths,n_materials,mean_tension=1):

    #Index interfaces :
    T = np.array(sorted(list(dict_lengths.keys())))
    kt = np.amax(T)+1
    Keys_t = T[:,0]*kt + T[:,1]
    reverse_map_t = dict(zip(Keys_t, np.arange(len(Keys_t))))
    
    #Index angles : 
    A = np.array(list(dict_rad.keys()))
    A = -np.sort(-A,axis=1)
    ka = np.amax(A)
    Keys_a = A[:,0]*ka**2 + A[:,1]*ka + A[:,2]
    Keys_a,index = np.unique(Keys_a,return_index=True)
    A=A[index]
    
    #MX = B
    #x : structure : [0,nm[ : tensions
    nm = len(T)
    nc = n_materials-1
    nj = len(A)
    size = (3*nj + 1)*(nm)
    B = np.zeros(3*nj + 1)
    B[-1]=nm*mean_tension  #Sum of tensions = number of membrane <=> mean of tensions = 1
    M = np.zeros((3*nj + 1, nm))
    
    
    #CLASSICAL 
    #YOUNG-DUPRÉ
    for i,key in enumerate(Keys_a):
        c = key//ka**2
        b = (key-(c*(ka**2)))//ka
        a = key%ka

        O_b = dict_rad[(a,b,c)] #O_b = abc angle
        O_a = dict_rad[(b,a,c)] #O_a = cab angle
        O_c = dict_rad[(a,c,b)]
        idx_tab = reverse_map_t[min(a,b)*kt+max(a,b)]
        idx_tbc = reverse_map_t[min(b,c)*kt+max(b,c)]
        idx_tac = reverse_map_t[min(a,c)*kt+max(a,c)]

        z1,z2,z3 = compute_z(O_a,O_b,O_c)
        f1,f2,f3,s1,s2,s3 = factors_z(z1,z2,z3)

        M[3*i,idx_tab]=-s3-f3    
        M[3*i,idx_tbc]=s3        
        M[3*i,idx_tac]=s3        

        M[3*i+1,idx_tab]=s1      
        M[3*i+1,idx_tbc]=-s1-f1  
        M[3*i+1,idx_tac]=s1      

        M[3*i+2,idx_tab]=s2      
        M[3*i+2,idx_tbc]=s2      
        M[3*i+2,idx_tac]=-s2-f2  

    M[3*nj,:nm]=1 #Enforces mean of tensions = 1 (see B[-1])

    return(M,B)


def build_matrix_tension_inv_cotan(dict_rad,dict_lengths,n_materials,mean_tension=1):

    #Index interfaces :
    T = np.array(sorted(list(dict_lengths.keys())))
    kt = np.amax(T)+1
    Keys_t = T[:,0]*kt + T[:,1]
    reverse_map_t = dict(zip(Keys_t, np.arange(len(Keys_t))))
    
    #Index angles : 
    A = np.array(list(dict_rad.keys()))
    A = -np.sort(-A,axis=1)
    ka = np.amax(A)
    Keys_a = A[:,0]*ka**2 + A[:,1]*ka + A[:,2]
    Keys_a,index = np.unique(Keys_a,return_index=True)
    A=A[index]
    
    #MX = B
    #x : structure : [0,nm[ : tensions
    nm = len(T)
    nc = n_materials-1
    nj = len(A)
    size = (3*nj + 1)*(nm)
    B = np.zeros(3*nj + 1)
    B[-1]=nm*mean_tension  #Sum of tensions = number of membrane <=> mean of tensions = 1
    M = np.zeros((3*nj + 1, nm))
    
    
    #CLASSICAL 
    #YOUNG-DUPRÉ
    for i,key in enumerate(Keys_a):
        c = key//ka**2
        b = (key-(c*(ka**2)))//ka
        a = key%ka

        O_b = dict_rad[(a,b,c)] #O_b = abc angle
        O_a = dict_rad[(b,a,c)] #O_a = cab angle
        O_c = dict_rad[(a,c,b)]
        idx_tab = reverse_map_t[min(a,b)*kt+max(a,b)]
        idx_tbc = reverse_map_t[min(b,c)*kt+max(b,c)]
        idx_tac = reverse_map_t[min(a,c)*kt+max(a,c)]

        z1,z2,z3 = compute_z(O_a,O_b,O_c)
        f1,f2,f3,s1,s2,s3 = factors_z(z1,z2,z3)

        M[3*i,idx_tab]=-1/s3-1/f3    
        M[3*i,idx_tbc]=1/f3        
        M[3*i,idx_tac]=1/f3        

        M[3*i+1,idx_tab]=1/f1      
        M[3*i+1,idx_tbc]=-1/s1-1/f1  
        M[3*i+1,idx_tac]=1/f1      

        M[3*i+2,idx_tab]=1/f2      
        M[3*i+2,idx_tbc]=1/f2      
        M[3*i+2,idx_tac]=-1/f2-1/s2  

    M[3*nj,:nm]=1 #Enforces mean of tensions = 1 (see B[-1])

    return(M,B)























def infer_tension_projection_yd(Mesh,mean_tension=1,check_topo_changes=False):
    dict_rad = Mesh.compute_angles()
    dict_lengths = Mesh.compute_lengths()
    n_materials = len(Mesh.faces)
    M,B = build_matrix_tension_projection_yd(dict_rad,dict_lengths,n_materials,mean_tension)
    x, resid, rank, sigma = linalg.lstsq(M, B)
    dict_tensions = {}
    T = np.array(sorted(list(dict_lengths.keys())))
    nm = len(T)    
    for i in range(nm): 
        dict_tensions[tuple(T[i])]=x[i]
        
    if check_topo_changes : 
        check_trijunctions(dict_rad,dict_lengths,dict_tensions)

    return(x,dict_tensions,resid)




def infer_tension_symmetrical_yd(Mesh,mean_tension=1,check_topo_changes=False):
    dict_rad = Mesh.compute_angles()
    dict_lengths = Mesh.compute_lengths()
    n_materials = len(Mesh.faces)
    M,B = build_matrix_tension_symmetrical_yd(dict_rad,dict_lengths,n_materials,mean_tension)
    x, resid, rank, sigma = linalg.lstsq(M, B)
    dict_tensions = {}
    T = np.array(sorted(list(dict_lengths.keys())))
    nm = len(T)    
    for i in range(nm): 
        dict_tensions[tuple(T[i])]=x[i]
        
    if check_topo_changes : 
        check_trijunctions(dict_rad,dict_lengths,dict_tensions)

    return(x,dict_tensions,resid)

def infer_tension_cotan(Mesh,mean_tension=1,check_topo_changes=False):
    dict_rad = Mesh.compute_angles()
    dict_lengths = Mesh.compute_lengths()
    n_materials = len(Mesh.faces)
    M,B = build_matrix_tension_cotan(dict_rad,dict_lengths,n_materials,mean_tension)
    x, resid, rank, sigma = linalg.lstsq(M, B)
    dict_tensions = {}
    T = np.array(sorted(list(dict_lengths.keys())))
    nm = len(T)    
    for i in range(nm): 
        dict_tensions[tuple(T[i])]=x[i]
        
    if check_topo_changes : 
        check_trijunctions(dict_rad,dict_lengths,dict_tensions)

    return(x,dict_tensions,resid)




def infer_tension_inv_cotan(Mesh,mean_tension=1,check_topo_changes=False):
    dict_rad = Mesh.compute_angles()
    dict_lengths = Mesh.compute_lengths()
    n_materials = len(Mesh.faces)
    M,B = build_matrix_tension_inv_cotan(dict_rad,dict_lengths,n_materials,mean_tension)
    x, resid, rank, sigma = linalg.lstsq(M, B)
    dict_tensions = {}
    T = np.array(sorted(list(dict_lengths.keys())))
    nm = len(T)    
    for i in range(nm): 
        dict_tensions[tuple(T[i])]=x[i]
        
    if check_topo_changes : 
        check_trijunctions(dict_rad,dict_lengths,dict_tensions)

    return(x,dict_tensions,resid)




def infer_tension_lamy(Mesh,mean_tension=1,check_topo_changes=False):
    dict_rad = Mesh.compute_angles()
    dict_lengths = Mesh.compute_lengths()
    n_materials = len(Mesh.faces)
    M,B = build_matrix_tension_lamy(dict_rad,dict_lengths,n_materials,mean_tension)
    x, resid, rank, sigma = linalg.lstsq(M, B)
    dict_tensions = {}
    T = np.array(sorted(list(dict_lengths.keys())))
    nm = len(T)    
    for i in range(nm): 
        dict_tensions[tuple(T[i])]=x[i]
        
    if check_topo_changes : 
        check_trijunctions(dict_rad,dict_lengths,dict_tensions)

    return(x,dict_tensions,resid)



def infer_tension_inv_lamy(Mesh,mean_tension=1,check_topo_changes=False):
    dict_rad = Mesh.compute_angles()
    dict_lengths = Mesh.compute_lengths()
    n_materials = len(Mesh.faces)
    M,B = build_matrix_tension_inv_lamy(dict_rad,dict_lengths,n_materials,mean_tension)
    x, resid, rank, sigma = linalg.lstsq(M, B)
    dict_tensions = {}
    T = np.array(sorted(list(dict_lengths.keys())))
    nm = len(T)    
    for i in range(nm): 
        dict_tensions[tuple(T[i])]=x[i]
        
    if check_topo_changes : 
        check_trijunctions(dict_rad,dict_lengths,dict_tensions)

    return(x,dict_tensions,resid)



def infer_tension_lamy_log(Mesh,mean_tension=1,check_topo_changes=False):
    dict_rad = Mesh.compute_angles()
    dict_lengths = Mesh.compute_lengths()
    n_materials = len(Mesh.faces)
    M,B = build_matrix_tension_lamy_log(dict_rad,dict_lengths,n_materials,mean_tension)
    x, resid, rank, sigma = linalg.lstsq(M, B)
    dict_tensions = {}
    T = np.array(sorted(list(dict_lengths.keys())))
    nm = len(T)    
    for i in range(nm): 
        dict_tensions[tuple(T[i])]=np.exp(x[i])
        
    if check_topo_changes : 
        check_trijunctions(dict_rad,dict_lengths,dict_tensions)

    return(x,dict_tensions,resid)



def plot_tension_inference(Mesh,dict_tensions=None,original_image = None): 
    if dict_tensions == None: 
        _,dict_tensions,_ = infer_tension(Mesh)
    Faces_list = Mesh.faces
    values = np.array(list(dict_tensions.values()))
    minv = np.amin(values)
    maxv = np.amax(values)
    dict_colors = {}
    for key in dict_tensions : 
        dict_colors[key] = plt.cm.jet((dict_tensions[key]-minv)/(maxv-minv))

    fig,ax = plt.subplots()

    Points_interface={}
    for f in Faces_list : 
        if not f.closed() : 
            continue

        he,success = f.find_triple_point_edge()
        if not success : 
            he = f.outer_component

        b=he

        while True : 
            Faces_key_unordered = (b.incident_face.attached['key'],b.twin.incident_face.attached['key'])
            Faces_key = (min(Faces_key_unordered),max(Faces_key_unordered))
            
            plt.plot([b.origin.y,b.destination.y],[b.origin.x,b.destination.x],color = dict_colors[Faces_key])
            b=b.next
            if (b.attached['key']==he.attached['key']) : 
                break
                
    ax.set_aspect('equal')
    plt.title("Surface tension inference")
    plt.gca().invert_yaxis()
    plt.xticks([])
    plt.yticks([])
    sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=minv, vmax=maxv))
    cbar = plt.colorbar(sm)
    cbar.set_label('Surface tensions')
    try: 
        plt.imshow(original_image,plt.cm.gray)
    except: 
        pass