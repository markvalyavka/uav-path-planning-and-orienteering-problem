from casadi import *
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def read_samples(filename):
    samples = None
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        header = next(csvreader)
        print(header)
        for row in csvreader:
            cols = []
            for c in row:
                #print(c)
                cols.append(float(c))
            if samples is None:
                samples = np.array(cols)
            else:
                samples = np.vstack((samples,np.array(cols)))
            #samples.append(cols)
    header_dic =  {x:int(i) for i,x in enumerate(header)}
    return samples, header_dic


max_jerk = np.array([[440,440,440]])

if __name__ == '__main__':
    N_points = 40
    ref,ref_dic = read_samples('../samples_equidistant_08.csv');

    opti = Opti()
    t = opti.variable(N_points-1)
    p = opti.variable(N_points,3)
    v = opti.variable(N_points,3)
    a = opti.variable(N_points,3)

    print(np.transpose(ref[0,ref_dic['p_x']:ref_dic['p_z']+1]))
    print(p[0,:])
    opti.subject_to(p[0,0] == ref[0,ref_dic['p_x']])
    opti.subject_to(p[0,1] == ref[0,ref_dic['p_y']])
    opti.subject_to(p[0,2] == ref[0,ref_dic['p_z']])
    opti.subject_to(v[0,0] == ref[0,ref_dic['v_x']])
    opti.subject_to(v[0,1] == ref[0,ref_dic['v_y']])
    opti.subject_to(v[0,2] == ref[0,ref_dic['v_z']])
    opti.subject_to(p[-1,0] == ref[N_points-1,ref_dic['p_x']])
    opti.subject_to(p[-1,1] == ref[N_points-1,ref_dic['p_y']])
    opti.subject_to(p[-1,2] == ref[N_points-1,ref_dic['p_z']])
    

    f = 0;
    for i in range(0,N_points-1):
        f+=t[i]
    G = []
    for i in range(0,N_points):
        opti.subject_to( a[i,1]*a[i,1]+a[i,1]*a[i,1]+a[i,1]*a[i,1] <= 32.0*32.0 )
    
    for i in range(1,N_points):
        t_part = t[i-1,0];
        opti.subject_to(p[i,:]==p[i-1,:]+v[i-1,:]*t_part+0.5*a[i-1,:]*t_part*t_part)
        opti.subject_to(v[i,:]==v[i-1,:]+a[i-1,:]*t_part)
        opti.subject_to(a[i,:]<=a[i-1,:]+max_jerk*t_part)
        opti.subject_to(a[i,:]>=a[i-1,:]-max_jerk*t_part)

        origt = ref[i,ref_dic['t']]-ref[i-1,ref_dic['t']]
        opti.subject_to(t_part>=0.8*origt)
        opti.subject_to(t_part<=1.0*origt)
        
        # probmin.Constraints.(const_name_p) = p(i,:)==p(i-1,:)+v(i-1,:)*t_part+0.5*a(i-1,:)*t_part^2;
        # probmin.Constraints.(const_name_v) = v(i,:)==v(i-1,:)+a(i-1,:)*t_part;
        # probmin.Constraints.(const_name_ajplus) = a(i,:)<=a(i-1,:)+max_jerk*t_part;
        # probmin.Constraints.(const_name_ajminus) = a(i,:)>=a(i-1,:)-max_jerk*t_part;
        
        
    opti.solver("ipopt")
    sol = opti.solve()
    print(sol.stats()["iter_count"])

    print("t is",sum(sol.value(t)))
    print("original time is ",ref[N_points-1,ref_dic['t']])

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(ref[0:N_points,ref_dic['p_x']], ref[0:N_points,ref_dic['p_y']], ref[0:N_points,ref_dic['p_z']], 'gray')
    plt.show()