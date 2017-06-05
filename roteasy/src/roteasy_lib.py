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

def rotate(cord,angles=(),axes='',unpacked=False,unpack=False,reference='r',xoff=0,yoff=0,zoff=0):
    """Rotate the point with coordinate 'cord' with an anticlokwise angle 'angles' around the axes 'axes'.
    The function return the new rotated coordinates of the points.

    :param cord: np.array (Nx3) or list with (x,y,z)
    :param angles: List of (anticlockwise) angle of rotation.
    :param axes: List of axes of rotation.
    :param unpacked: If True, coord are on the form (x,y,z), if False cord are a Nx3 array with x in the 0-col, y 1-col and z 2-col
    :param unpack: If True return the new cord in the form (xn,yn,zn), if False return the coord in a Nx3 array with x in the 0-col, y 1-col and z 2-col
    :param reference: rh if the system is right-handed, lh if the system is left-handed.
    :param xoff: X-offset after the rotation
    :param yoff: Y-offset after the rotation
    :param zoff: Z-offset after the rotation
    :return:
    """

    if reference[0].lower()== 'r': cost2 = 1
    elif reference[0].lower()== 'l':cost2 = -1
    else: raise ValueError('Wrong frame formalism')



    angles=np.array(angles)*cost2


    cord=np.array(cord)

    if unpacked: xo, yo, zo = cord
    else: xo,yo,zo = cord.T



    R=rotation_matrix(angles=angles,axes=axes)



    xn = R[0,0]*xo + R[1,0]*yo + R[2,0]*zo - xoff
    yn = R[0,1]*xo + R[1,1]*yo + R[2,1]*zo - yoff
    zn = R[0,2]*xo + R[1,2]*yo + R[2,2]*zo - zoff



    ret = (xn,yn,zn)

    if unpack: return ret
    else: return np.array(ret).T

def rotate_frame(cord,angles=(),axes='',unpacked=False,unpack=False,reference='r',xoff=0,yoff=0,zoff=0):
    """
    Rotate the frame of reference with an anticlokwise angle 'angles' around the axes 'axes'.
    The function return the coordinates of the points 'cord' in the new frame of reference.
    The frame is rotated with angles[0] around the axes[0], then it is rotated of the angles[1] around the new axes axes[1] and so on...
    :param cord: np.array (Nx3) or list with (x,y,z)
    :param angles: List of (anticlockwise) angle of rotation.
    :param axes: List of axes of rotation.
    :param unpacked: If True, coord are on the form (x,y,z), if False cord are a Nx3 array with x in the 0-col, y 1-col and z 2-col
    :param unpack: If True return the new cord in the form (xn,yn,zn), if False return the coord in a Nx3 array with x in the 0-col, y 1-col and z 2-col
    :param reference: rh if the system is right-handed, lh if the system is left-handed.
    :param xoff: X-offset after the rotation
    :param yoff: Y-offset after the rotation
    :param zoff: Z-offset after the rotation
    :return:
    """


    angles=-1*np.array(angles)


    return  rotate(cord=cord,angles=angles,axes=axes,unpacked=unpacked,unpack=unpack,reference=reference,xoff=xoff,yoff=yoff,zoff=zoff)

def align_frame(cord,pos_vec,ax='x',cartesian=True,spherical=False,cylindrical=False,unpacked=False,unpack=False,reference='r',xoff=0,yoff=0,zoff=0):
    """
    Rotate the frame of reference to align the axis specied in ax to the vector specied in pos_vec.
    The function return the position of the particles with cord 'cord' in the new frame of reference.
    :param cord: np.array (Nx3) or list with (x,y,z)
    :param pos_vec:  Vector containing the coordinate of the reference vector (see below)
    :param ax: ax to align to Vector
    :param cartesian: the coordinates of pos_vec are cartesian (x,y,z)
    :param spherical: the coordinates of pos_vec are spherical  (phi (azimuth), theta(declination)) in degree
    :param cylindrical:  the coordinates of pos_vec are
    :param xoff: X-offset after the rotation
    :param yoff: Y-offset after the rotation
    :param zoff: Z-offset after the rotation
    :return:
    """

    if cartesian: l,b=cartesian_to_spherical(pos_vec)
    elif spherical: l,b=pos_vec
    elif cylindrical: l,b=cylindrical_to_spherical(pos_vec)
    else: ValueError('At least one must be true among keyword cartesian,spherical and cylindrical')

    if reference[0].lower()== 'r': cost2 = 1
    elif reference[0].lower()== 'l':cost2 = -1
    else: raise ValueError('Wrong frame formalism')

    if ax=='x':
        angles=(l*cost2,-b*cost2)
        axes='zy'
    elif ax=='-x':
        angles=(l*cost2, (-b+180)*cost2 )
        axes='zy'
    elif ax=='y':
        angles=( (270+l)*cost2, b*cost2)
        axes='zx'
    elif ax=='-y':
        angles=( (270+l)*cost2, (b+180)*cost2)
        axes='zx'
    elif ax=='z':
        angles=(l*cost2, (-b+90)*cost2)
        axes='zy'
    elif ax=='-z':
        angles=(l*cost2, (-b+90+180)*cost2)
        axes='zy'
    else:
        raise ValueError('ax %s unknown'%ax)



    return rotate_frame(cord=cord,angles=angles,axes=axes,unpacked=unpacked,unpack=unpack,reference=reference,xoff=xoff,yoff=yoff,zoff=zoff)

def rotate2D(cord,angle,unpacked=False,unpack=False,reference='r',xoff=0,yoff=0):
    """
    Rotate a point in a X-Y frame of reference. The rotation is around the Z axis.
    :param cord: np.array (Nx2) or list with (x,y)
    :param angles: (anticlockwise) angle of rotation.
    :param unpacked: If True, coord are on the form (x,y), if False cord are a Nx3 array with x in the 0-col, y 1-col
    :param unpack: If True return the new cord in the form (xn,yn), if False return the coord in a Nx2 array with x in the 0-col, y 1-col
    :param reference: rh if the system is right-handed, lh if the system is left-handed.
    :param xoff: X-offset after the rotation
    :param yoff: Y-offset after the rotation
    :param zoff: Z-offset after the rotation
    :return:
    """
    angles=(angle,)
    ax='z'
    zoff=0


    if unpacked: xo, yo = cord
    else: xo,yo = cord.T

    zo=np.zeros_like(xo,dtype=float)

    xn,yn,_=rotate(cord=(xo,yo,zo),angles=angles,axes=ax,unpacked=True,unpack=True,reference=reference,xoff=xoff,yoff=yoff,zoff=zoff)

    if unpack: return (xn,yn)
    else: return np.array((xn,yn)).T

def rotate_frame2D(cord,angle,unpacked=False,unpack=False,reference='r',xoff=0,yoff=0):
    """
    Rotate  a X-Y frame of reference. The rotation is around the Z axis.
    The function return the coordinates of the points 'cord' in the new frame of reference.
    :param cord: np.array (Nx2) or list with (x,y)
    :param angles: (anticlockwise) angle of rotation.
    :param unpacked: If True, coord are on the form (x,y), if False cord are a Nx3 array with x in the 0-col, y 1-col
    :param unpack: If True return the new cord in the form (xn,yn), if False return the coord in a Nx2 array with x in the 0-col, y 1-col
    :param reference: rh if the system is right-handed, lh if the system is left-handed.
    :param xoff: X-offset after the rotation
    :param yoff: Y-offset after the rotation
    :param zoff: Z-offset after the rotation
    :return:
    """
    return rotate2D(cord,-angle,unpacked=unpacked,unpack=unpack,reference=reference,xoff=xoff,yoff=yoff)

def align_frame2D(cord,pos_vec,ax='x',cartesian=True, polar=False,unpacked=False,unpack=False,reference='r',xoff=0,yoff=0,zoff=0):
    """
    Rotate a 2D X-Y frame of reference to align the axis specied in ax to the vector specied in pos_vec.
    The function return the position of the particles with cord 'cord' in the new frame of reference.
    :param cord: np.array (Nx2) or list with (x,y)
    :param pos_vec:  Vector containing the coordinate of the reference vector (see below)
    :param ax: ax to align to Vector (x, or y)
    :param cartesian: the coordinates of pos_vec are cartesian (x,y)
    :param polar: the aximuthal angle in degree
    :param xoff: X-offset after the rotation
    :param yoff: Y-offset after the rotation
    :return:
    """

    if cartesian: angle,_=cartesian_to_spherical((pos_vec[0],pos_vec[1],0))
    elif polar: angle=pos_vec
    else: ValueError('At least one must be true among keyword cartesian and polar')

    if reference[0].lower()== 'r': cost2 = 1
    elif reference[0].lower()== 'l':cost2 = -1
    else: raise ValueError('Wrong frame formalism')

    if ax=='x': angle=angle*cost2
    elif ax=='-x': angle=(angle+180)*cost2
    elif ax=='y': angle=(270+angle)*cost2
    elif ax=='-y': angle=(270+angle+180)*cost2
    else: raise ValueError('ax %s unknown'%ax)

    return  rotate_frame2D(cord=cord,angle=angle,unpacked=unpacked,unpack=unpack,reference=reference,xoff=xoff,yoff=yoff)

def rotate_euler(cord,alpha=0,beta=0,gamma=0,axes='zxz',unpacked=False,unpack=False,reference='r',xoff=0,yoff=0,zoff=0):
    """
    Rotate the frame of reference using the Euler convention: a first rotation of alpha round the original axes[0], then
    a rotation of angle beta around the new axes[1] and finally a rotation of gamma around the final axis axes[2].
    :param cord: np.array (Nx3) or list with (x,y,z)
    :param alpha: The first anticlokwise rotation angle
    :param beta: The second anticlokwise rotation angle
    :param gamma: The third anticlokwise rotation angle
    :param axes: axes of rotation
    :param unpacked: If True, coord are on the form (x,y,z), if False cord are a Nx3 array with x in the 0-col, y 1-col and z 2-col
    :param unpack: If True return the new cord in the form (xn,yn,zn), if False return the coord in a Nx3 array with x in the 0-col, y 1-col and z 2-col
    :param reference: rh if the system is right-handed, lh if the system is left-handed.
    :param xoff: X-offset after the rotation
    :param yoff: Y-offset after the rotation
    :param zoff: Z-offset after the rotation
    :return:
    """
    angles=(alpha,beta,gamma)

    return rotate_frame(cord=cord,angles=angles,axes=axes,unpack=unpack,unpacked=unpacked,reference=reference,xoff=xoff,yoff=yoff,zoff=zoff)


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

