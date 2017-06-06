from __future__ import division, print_function
import numpy as np
from math import cos,sin, sqrt, atan2, asin

angle_to_rad=np.pi/180.
rad_to_angle=1/angle_to_rad

def cartesian_to_spherical(pos):
    x,y,z=pos
    d=sqrt(x*x+y*y+z*z)

    l=atan2(y,x)*rad_to_angle

    b=asin(z/d)*rad_to_angle

    return l,b

def cylindrical_to_spherical(pos):

    R,z,l=pos

    b = asin(z / sqrt(R*R+z*z))*rad_to_angle

    return l,b

def rotation_matrix_x(angle):
    """Calculate the matrix that defines the anticlockwise rotation of a point around the x axis:

    .. math::
        :nowrap:

            \\begin{bmatrix}
                1       & 0 &  0 \\\\
                0       & \\cos \\theta  & \\sin \\theta  \\\\
                0      & -\\sin \\theta & \\cos \\theta
            \\end{bmatrix}

    :param angle: (anticlockwise) rotation angle ( :math:`\\theta`) in degree
    :return: a 3x3 numpy array containing the rotation matrix (see above).
    """

    angle_rad=angle*angle_to_rad

    ca=cos(angle_rad)
    sa=sin(angle_rad)

    return np.array([[1, 0, 0], [0, ca, sa], [0, -sa, ca]])

def rotation_matrix_y(angle):
    """Calculate the matrix that defines the anticlockwise rotation of a point around the y axis:

    .. math::
        :nowrap:

            \\begin{bmatrix}
                \\cos \\theta        & 0 &  -\\sin \\theta \\\\
                0       & 1  & 0\\\\
                \\sin \\theta      & 0 & \\cos \\theta
            \\end{bmatrix}

    :param angle: (anticlockwise) rotation angle ( :math:`\\theta`) in degree
    :return: a 3x3 numpy array containing the rotation matrix (see above).
    """

    angle_rad=angle*angle_to_rad

    cb=cos(angle_rad)
    sb=sin(angle_rad)

    return np.array([[cb, 0, -sb], [0, 1, 0], [sb, 0, cb]])

def rotation_matrix_z(angle):
    """Calculate the matrix that defines the anticlockwise rotation of a point around the z axis:

    .. math::
        :nowrap:

            \\begin{bmatrix}
                \\cos \\theta        & \\sin \\theta & 0 \\\\
                -\\sin \\theta        & \\cos \\theta  & 0\\\\
                0      & 0 & 1
            \\end{bmatrix}

    :param angle: (anticlockwise) rotation angle ( :math:`\\theta`) in degree
    :return: a 3x3 numpy array containing the rotation matrix (see above).
    """

    angle_rad=angle*angle_to_rad

    cg=cos(angle_rad)
    sg=sin(angle_rad)

    return  np.array([[cg,sg,0],[-sg,cg,0],[0,0,1]])

dic_rot={'x':rotation_matrix_x,'y':rotation_matrix_y,'z':rotation_matrix_z}

def rotation_matrix(angles=(),axes=''):
    """Core function to obtain the final rotation matrix :math:`R`.

    .. math::

        R(\\theta_1,\\theta_2,....,\\theta_N)=R_1(\\theta_1) \cdot R_2(\\theta_2) \cdot ..... \cdot R_N(\\theta_N)


    where :math:`\\theta_i` is the anticlockwise rotation angle around the axis :math:`x_i`
    e.g.::

        angles=(30,60,80)
        axes='xyz'

        R=rotation_matrix(angles=angles,axes=axes)

    In this case R is the rotation matrix representing a rotation of 30 degree around the
    x axis, then a rotation of 60 degree around the new y axis and finally a rotation of 80 degree
    around the final z axis.

    :param angles: (anticlockwise) rotation angles (list, tuple or numpy array).
    :param axes: str representing the rotation axes.
    :return: final rotation matrix (see above).
    """

    nangles=len(angles)
    naxes=len(axes)

    if (nangles==0) or (naxes==0): raise ValueError('Number of angles and/or rotation axes can not be 0 in rotation_matrix of module roteasy')
    elif nangles>naxes: raise ValueError('More angles than axes of rotation specified in rotation_matrix of module roteasy')
    elif naxes>nangles: print('Warning: More axes than angles, exceding axes skipped')

    R=np.array([[1.,0,0],[0,1.,0],[0,0,1.]])

    for i in range(nangles):
        angle=angles[i]
        axis=axes[i]
        rot_mat=dic_rot[axis](angle)

        R=np.dot(R,rot_mat)

    return R




"""
#TEST
#row convention, the following functions perform the same rotation but following
#the row convention (x,y,z).T= R (x.y,z).T NB both the singl matrix and the order of the multiplicaiton change
#They are not used!

def rotation_matrix_x(angle):

    angle_rad=angle*angle_to_rad

    ca=cos(angle_rad)
    sa=sin(angle_rad)

    return np.array([[1, 0, 0], [0, ca, -sa], [0, sa, ca]])

def rotation_matrix_y(angle):

    angle_rad=angle*angle_to_rad

    cb=cos(angle_rad)
    sb=sin(angle_rad)

    return np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]])

def rotation_matrix_z(angle):

    angle_rad=angle*angle_to_rad

    cg=cos(angle_rad)
    sg=sin(angle_rad)

    return  np.array([[cg,-sg,0],[sg,cg,0],[0,0,1]])

dic_rot={'x':rotation_matrix_x,'y':rotation_matrix_y,'z':rotation_matrix_z}

def rotation_matrix(angles=(),axes=''):


    nangles=len(angles)
    naxes=len(axes)

    if (nangles==0) or (naxes==0): raise ValueError('Number of angles and/or rotation axes can not be 0 in rotation_matrix of module roteasy')
    elif nangles>naxes: raise ValueError('More angles than axes of rotation specified in rotation_matrix of module roteasy')
    elif naxes>nangles: print('Warning: More axes than angles, exceding axes skipped')

    R=np.array([[1.,0,0],[0,1.,0],[0,0,1.]])

    for i in range(nangles-1,-1,-1):
        angle=angles[i]
        axis=axes[i]
        rot_mat=dic_rot[axis](angle)

        R=np.dot(rot_mat,R)

    return R

def rotate(cord,angles=(),axes='',unpacked=False,unpack=False,reference='r',xoff=0,yoff=0,zoff=0):


    if reference[0].lower()== 'r': cost2 = 1
    elif reference[0].lower()== 'l':cost2 = -1
    else: raise ValueError('Wrong frame formalism')



    angles=np.array(angles)*cost2


    cord=np.array(cord)

    if unpacked: xo, yo, zo = cord
    else: xo,yo,zo = cord.T



    R=rotation_matrix(angles=angles,axes=axes)



    xn = R[0,0]*xo + R[0,1]*yo + R[0,2]*zo - xoff
    yn = R[1,0]*xo + R[1,1]*yo + R[1,2]*zo - yoff
    zn = R[2,0]*xo + R[2,1]*yo + R[2,2]*zo - zoff



    ret = (xn,yn,zn)

    if unpack: return ret
    else: return np.array(ret).T
"""

