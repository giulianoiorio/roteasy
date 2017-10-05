from .src.roteasy_lib import *


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


    return rotate_core(cord=cord,angles=angles,axes=axes,unpacked=unpacked,unpack=unpack,reference=reference,xoff=xoff,yoff=yoff,zoff=zoff,change_reference=False)

def rotate_frame(cord,angles=(),axes='',unpacked=False,unpack=False,reference='r',xoff=0,yoff=0,zoff=0,change_reference=False):
    """Rotate the frame of reference with an anticlokwise angle 'angles' around the axes 'axes'.
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
    :param change_reference: If 'x' 'y' or 'z' change the final frame of reference (after all the rotations) to the convention opposite to the one indicating in reference. If '' or False, nothing change.
    :return:
    """


    angles=-1*np.array(angles)


    return rotate_core(cord=cord,angles=angles,axes=axes,unpacked=unpacked,unpack=unpack,reference=reference,xoff=xoff,yoff=yoff,zoff=zoff,change_reference=change_reference)

def align_frame(cord,pos_vec,ax='x',cartesian=True,spherical=False,cylindrical=False,unpacked=False,unpack=False,reference='r',xoff=0,yoff=0,zoff=0,change_reference=False):
    """Rotate the frame of reference to align the axis specied in ax to the vector specied in pos_vec.
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
    :param change_reference: If 'x' 'y' or 'z' change the final frame of reference (after all the rotations) to the convention opposite to the one indicating in reference. If '' or False, nothing change.
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
        angles=(l*cost2, (-b+180)*cost2)
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



    return rotate_frame(cord=cord,angles=angles,axes=axes,unpacked=unpacked,unpack=unpack,reference=reference,xoff=xoff,yoff=yoff,zoff=zoff,change_reference=change_reference)

def rotate2D(cord,angle,unpacked=False,unpack=False,reference='r',xoff=0,yoff=0):
    """Rotate a point in a X-Y frame of reference. The rotation is around the Z axis.

    :param cord: np.array (Nx2) or list with (x,y)
    :param angles: (anticlockwise) angle of rotation.
    :param unpacked: If True, coord are on the form (x,y), if False cord are a Nx3 array with x in the 0-col, y 1-col
    :param unpack: If True return the new cord in the form (xn,yn), if False return the coord in a Nx2 array with x in the 0-col, y 1-col
    :param reference: rh if the system is right-handed, lh if the system is left-handed.
    :param xoff: X-offset after the rotation
    :param yoff: Y-offset after the rotation
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

def rotate_frame2D(cord,angle,unpacked=False,unpack=False,reference='r',xoff=0,yoff=0,change_reference=False):
    """Rotate  a X-Y frame of reference. The rotation is around the Z axis.
    The function return the coordinates of the points 'cord' in the new frame of reference.

    :param cord: np.array (Nx2) or list with (x,y)
    :param angles: (anticlockwise) angle of rotation.
    :param unpacked: If True, coord are on the form (x,y), if False cord are a Nx3 array with x in the 0-col, y 1-col
    :param unpack: If True return the new cord in the form (xn,yn), if False return the coord in a Nx2 array with x in the 0-col, y 1-col
    :param reference: rh if the system is right-handed, lh if the system is left-handed.
    :param xoff: X-offset after the rotation
    :param yoff: Y-offset after the rotation
    :param change_reference: If 'x' 'y' or 'z' change the final frame of reference (after all the rotations) to the convention opposite to the one indicating in reference. If '' or False, nothing change.
    :return:
    """

    angles=(angle,)
    ax='z'
    zoff=0


    if unpacked: xo, yo = cord
    else: xo,yo = cord.T

    zo=np.zeros_like(xo,dtype=float)

    xn,yn,_=rotate(cord=(xo,yo,zo),angles=angles,axes=ax,unpacked=True,unpack=True,reference=reference,xoff=xoff,yoff=yoff,zoff=zoff,change_reference=change_reference)

    if unpack: return (xn,yn)
    else: return np.array((xn,yn)).T

def align_frame2D(cord,pos_vec,ax='x',cartesian=True, polar=False,unpacked=False,unpack=False,reference='r',xoff=0,yoff=0,change_reference=False):
    """Rotate a 2D X-Y frame of reference to align the axis specied in ax to the vector specied in pos_vec.
    The function return the position of the particles with cord 'cord' in the new frame of reference.

    :param cord: np.array (Nx2) or list with (x,y)
    :param pos_vec:  Vector containing the coordinate of the reference vector (see below)
    :param ax: ax to align to Vector (x, or y)
    :param cartesian: the coordinates of pos_vec are cartesian (x,y)
    :param polar: the aximuthal angle in degree
    :param xoff: X-offset after the rotation
    :param yoff: Y-offset after the rotation
    :param change_reference: If 'x' 'y' or 'z' change the final frame of reference (after all the rotations) to the convention opposite to the one indicating in reference. If '' or False, nothing change.
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

    return  rotate_frame2D(cord=cord,angle=angle,unpacked=unpacked,unpack=unpack,reference=reference,xoff=xoff,yoff=yoff,change_reference=change_reference)

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





