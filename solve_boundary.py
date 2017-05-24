#coding:utf-8

import numpy as np

def solve_boundary(a,xb,bf):
    xf=-bf
    bf_index=np.where(bf==True)[0]
    xf_index=np.where(xf==True)[0]
    b0=xb[bf_index]
    x0=xb[xf_index]

    a_x=a[xf][:,xf]
    a_xc_row=a[bf][:,xf]
    a_xc_column=a[xf][:,bf]
    a_b=a[bf][:,bf]

    b1=b0-np.dot(a_xc_row,x0)
    x1=np.linalg.solve(a_b,b1)
    b2=np.dot(a_xc_column,x1)+np.dot(a_x,x0)

    tx=[]
    tx2=None
    tb=[]
    tb2=None
    for i in range(len(bf)):
        if xf[i]==True:
            tx2,x0=x0[0],x0[1:]
            tb2,b2=b2[0],b2[1:]
        else:
            tx2,x1=x1[0],x1[1:]
            tb2,b0=b0[0],b0[1:]
        tx.append(tx2)
        tb.append(tb2)
    b=np.array(tb)
    x=np.array(tx)

    return b,x
    
if __name__=="__main__":

    """
    a=np.arange(1,17).reshape(4,4)
    x=np.array([4,3,2,1])
    b=np.array([20,60,100,140])
    """
    a=np.arange(1,17).reshape(4,4)
    xb=np.array([20,60,2,1])
    bf=np.array([True,True,False,False])

    b,x=solve_boundary(a,xb,bf)

