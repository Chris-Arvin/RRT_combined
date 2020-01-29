import numpy as np
from cvxopt import matrix, solvers
import pprint
import cvxopt
import tkinter as tk
from rrt_line import *
# from geometry_msgs import *

from scipy import optimize as op

'''
rrt_crooked___update log
1. in this py file, we can receive a linear path from rrt_line.py
and then we try to smooth the path

rrt_crooked1___update log
1. we figure out a bug about p:在solver中，输入的必须是矩阵，而矩阵的前身应该是数组，在ｐ这种ｎｘ１的形式中，不能用list进行append，否则会导致形式错误(shape不是ｎｘ１，而是n,)．
这里应该用np.ones，即直接在array上进行操作；注意，np.zeros生成的数组，里面全是ｉｎｔ型，如果想生成数组形，要用np.ones*0.0

rrt_crooked2___updatelog
1. add the function that when the springed crook meet with collisions, we add a new point between two old points. in this way, it can greatly escape from stacking collisions
'''

def advance():
    #find t_every and t_to
    global T
    T = 1
    global t_every
    t_every=[]
    Distance_all = 0
    for i in range(len(list) - 1):
        dis = np.sqrt((list[i + 1][1] - list[i][1]) ** 2 + (list[i + 1][0] - list[i][0]) ** 2)
        Distance_all += dis
    for i in range(len(list) - 1):
        dis = np.sqrt((list[i + 1][1] - list[i][1]) ** 2 + (list[i + 1][0] - list[i][0]) ** 2)
        t_every.append(dis / Distance_all * T)

    global t_to
    t_to=[]
    t_to.append(0)
    for i in range(len(t_every)):
        t_to.append(t_to[i]+t_every[i])

#define Q
def Q_temp(length=6):
    global T
    q_temp=np.ones([length,length])*0.0
    for i in range(length):
        for j in range(length):
            if i>=3 and j>=3:
                q_temp[i][j]=(i+1)*(i)*(i-1)*(i-2)*(j+1)*(j)*(j-1)*(j-2)/(i+j-5)*T**(i+j-5)
    return q_temp

def def_Q_all():
    global Q_all
    Q_all = []
    for i in range(k * (n + 1)):
        q_t = []
        for j in range(k * (n + 1)):
            q_t.append(0.0)
        Q_all.append(q_t)

    for i in range(k):
        q_temp = Q_temp()
        num_q = i * (n + 1)
        for j in range(n + 1):
            for u in range(n + 1):
                Q_all[num_q + j][num_q + u] = q_temp[j][u]
    Q_all = np.mat(Q_all)
    print('Q_all is done:')
    print(Q_all.shape)
    print('-'*60)



#define constrains
def p_constrain(list_p,a0=0,v0=0,ak=0,vk=0):
    global p,M
    #define constrains: firstly equal constrains
    #p
    p=np.ones([4*k+2,1])*0.0
    p[0][0]=list_p[0]
    p[1][0]=a0
    p[2][0]=v0
    # p.append([list_p[0]])
    # p.append([a0])
    # p.append([v0])

    for i in range(1,k):
        p[i+2][0]=list_p[i]

    p[k+2][0]=list_p[k]
    p[k+3][0]=ak
    p[k+4][0]=vk
    # p.append([list_p[k]])
    # p.append([ak])
    # p.append([vk])
    print('p is done')
    print('k+5的值:',k+5)
    print('p的总长度:',len(p),'4k+2的值',4*k+2)
    p=np.array(p)
    print(p.shape)
    print('-' * 60)

    M=np.ones([4*k+2,(n+1)*k])*0.0
    #M1
    M1=np.ones([k-1+3+3,(n+1)*k])*0.0
    #p0,v0,a0
    for i in range(n+1):
        M1[0][i]=t_to[0]**i
        if i==n:
            continue
        M1[1][i+1]=(i+1)*t_to[0]**i
        if i==n-1:
            continue
        M1[2][i+2]=(i+2)*(i+1)*t_to[0]**i

    #middle p
    for j in range(3,k+2):
        for i in range(n + 1):
            M1[j][i+(j-2)*(n+1)] = t_to[j-2] ** i

    #pk,vk,ak
    for i in range(n+1):
        M1[k+2][i+(k-1)*(n+1)]=t_to[-1]**i
        if i==n:
            continue
        M1[k+3][i+1+(k-1)*(n+1)]=(i+1)*t_to[-1]**i
        if i==n-1:
            continue
        M1[k+4][i+2+(k-1)*(n+1)]=(i+2)*(i+1)*t_to[-1]**i

    print('M1 is done')
    print(M1.shape)
    print('-' * 60)



    #definew constrains:secondly unequal constrains
    M2=np.ones([3*k-3,(n+1)*k])*0.0
    j=0
    l=0
    while l<3*(k-1):
        for i in range(n+1):
            M2[l][j*(n+1)+i]=t_to[j+1]**i
            M2[l][j * (n + 1) + i+(n+1)] = -t_to[j+1] ** i
            if i == n:
                continue
            M2[l+1][j*(n+1)+i + 1] = (i + 1) * t_to[j+1] ** i
            M2[l+1][j*(n+1)+i + 1+(n+1)] = -(i + 1) * t_to[j+1] ** i
            if i == n - 1:
                continue
            M2[l+2][j*(n+1)+i + 2] = (i + 2) * (i + 1) * t_to[j+1] ** i
            M2[l+2][j*(n+1)+i + 2+(n+1)] = - (i + 2) * (i + 1) * t_to[j+1] ** i

        j+=1
        l+=3
    print('M2 is done')
    print(M2.shape)
    print('-' * 60)

    #combine M1 and M2 to M
    #M=np.ones([4*k+2,(n+1)*k])*0.0
    for i in range(4*k+2):
        for j in range((n+1)*k):
            if i<k+5:
                M[i][j]=M1[i][j]
            else:
                M[i][j]=M2[i-k-5][j]

    # for i in range(26):
    #     print(M[i])
    #     if i<11:
    #         print(M1[i])
    #     else:
    #         print(M2[i-11])
    #     print('-'*60)
    print('M is done')
    print(M.shape)
    print('-' * 60)

    # figure_out()
    global Q_all
    q = np.ones([len(Q_all), 1])
    q = np.mat(q)
    # cvxopt_solve_qp(P=Q_all, q=q,A=M, b=p)
    # q=np.transpose(q)
    # Q_all=np.transpose(Q_all)
    # M=np.transpose(M)
    # p=np.transpose(p)
    q = matrix(q)
    Q_all = matrix(Q_all)
    M = matrix(M)
    p = matrix(p)
    # p=matrix([3.0,4.0,0.0,0.0,0.0,0.0,0,0,0,0])
    # print(Q_all.shape)
    # print(M.shape)
    # print(p.shape)
    result = solvers.qp(P=Q_all, q=q, A=M, b=p)
    return result['x']



def figure_out():

    advance()
    def_Q_all()
    lama_x=p_constrain(list_x)

    advance()
    def_Q_all()
    lama_y=p_constrain(list_y)
    global x,y
    x=[]
    y=[]
    global T
    global t_to,t_every


    #distrubute every segment into 20 points
    for t_a in range(len(t_to)-1):
        time_list = np.linspace(t_to[t_a], t_to[t_a+1], 20,endpoint=True)
        for j in time_list:
            m = np.ones([len(lama_x), 1]) * 0.0
            m[0+t_a*6]=1
            m[1+t_a*6]=j
            m[2+t_a*6] = j ** 2
            m[3+t_a*6] = j ** 3
            m[4+t_a*6] = j ** 4
            m[5+t_a*6] = j ** 5
            x.append(np.dot(np.transpose(m), lama_x)[0][0])
            y.append(np.dot(np.transpose(m), lama_y)[0][0])
            if rrt_agent.col_map[int(np.dot(np.transpose(m), lama_x)[0][0])][int(np.dot(np.transpose(m), lama_y)[0][0])]>50:
                l=list.copy()
                b1=int((list[t_a][0]+list[t_a+1][0])/2)
                b2=int((list[t_a][1]+list[t_a+1][1])/2)
                l.insert(t_a+1,[b1,b2])
                return False,l

    print('-'*60)
    print('x:')
    print(x)
    print('y:')
    print(y)
    print('-' * 60)
    return True,1



def draw():
    for yyy in range(1, rrt_agent.height - 1):
        for xxx in range(1, rrt_agent.width - 1):
            if rrt_agent.col_map[xxx][yyy] > 50:
                canvas.create_rectangle(xxx - 1, yyy - 1, xxx + 1,
                                              yyy + 1,
                                              fill='black')
    for i in range(len(rrt_agent.path_end)):
        canvas.create_rectangle(rrt_agent.path_end[i].col - 2, rrt_agent.path_end[i].row - 2, rrt_agent.path_end[i].col + 2,
                                      rrt_agent.path_end[i].row + 2,
                                      fill='green')

    for i in range(len(x)-1):
        print(x[i],y[i])
        canvas.create_line(int(x[i]),int(y[i]),int(x[i+1]),int(y[i+1]),fill='red')
        canvas.update()

if __name__ == '__main__':
    rrt_agent = rrt()
    rrt_agent.init_map()
    list=rrt_agent.path_xy
    n = 5
    while True:
        list_x = []
        list_y = []
        k = len(list) - 1
        for i in list:
            list_x.append(i[0])
            list_y.append(i[1])
        is_success,ll=figure_out()
        if is_success:
            break
        list=ll


    window = tk.Tk()
    window.title('rrt')
    window.geometry('%dx%d' % (800,900))
    canvas=tk.Canvas(window, bg='white', height=800, width=800)
    b = tk.Button(window, text='draw', command=draw)
    canvas.place(x=0,y=0,anchor='nw')
    b.place(x=400,y=850,anchor='nw')
    window.mainloop()

