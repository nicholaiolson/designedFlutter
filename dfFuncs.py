# -*- coding: utf-8 -*-
"""
Nicholai Olson
230214

Purpose: Function library for Designed Flutter project
"""
import numpy as np
from shapely import geometry
import drangSectionProps 
import matplotlib.pyplot as plt

def readAirfoil(file, headerLines): 
    """
    Function to read an airfoil from a standard .dat text file airfoil format. 

    Parameters
    ----------
    file : str
        Selig format airfoil .dat file format (Xfoil format)
    headerLines : int
        Number of lines containing header information for the .dat file

    Returns
    -------
    Airfoil data in 2D Numpy array

    """
    
    
    #Open text file
    airfoilFile = open(file, "r")
    fileContent = airfoilFile.readlines()
    
    #Loop through and collect airfoil points
    points = []
    for line in fileContent[headerLines:]:
        points.append(line.split())
    
    #Convert to numpy array
    airfoilArray = np.array(points,'float')
    
    return airfoilArray
    
def sectionPropsCheck(profile):
    """
    Computes seciton properties of a basic profile of uniform thickness, 
    such as an airfoil. 

    Parameters
    ----------
    profile : 2D np.array
        Numpy array of profile such as Selig format airfoil data
    Returns
    -------
    Section properties: 
        Ixx, Iyy, Ixy, A, x_bar, y_bar

    """
    
    #Create Shapely LineString objects for scaling and offsetting
    OML = geometry.polygon.LineString(profile)
    
    IXX_OML_arr = np.zeros([len(OML.xy[0]),1])
    IYY_OML_arr = np.zeros([len(OML.xy[0]),1])
    IXY_OML_arr = np.zeros([len(OML.xy[0]),1])
    x_bar_num_arr = np.zeros([len(OML.xy[0]),1])
    y_bar_num_arr = np.zeros([len(OML.xy[0]),1])
    dA = np.zeros([len(OML.xy[0]),1])
    for iPoint, point in enumerate(OML.xy[0][1:]): 
        
        #Moments of area of the filled shape about the origin
        #Perform Ixx area calc about Y axis for better discretization relative to X axis
        dA[iPoint,0] = (OML.xy[1][iPoint]-OML.xy[1][iPoint-1])*(OML.xy[0][iPoint]+OML.xy[0][iPoint-1])/2
        x_dA = (OML.xy[0][iPoint]+OML.xy[0][iPoint-1])/2/2
        y_dA = (OML.xy[1][iPoint]+OML.xy[1][iPoint-1])/2
        IXX_OML_arr[iPoint,0] = abs(y_dA**2*dA[iPoint,0])
        y_bar_num_arr[iPoint,0] = dA[iPoint,0]*y_dA
        
        #Perform Iyy area calc about X axis for better discritazation relative to Y axis
        #Also compute IXY about this axis as it'll have less numerical accuracy for an airfoil
        dA[iPoint,0] = (OML.xy[1][iPoint]+OML.xy[1][iPoint-1])/2*(OML.xy[0][iPoint]-OML.xy[0][iPoint-1])
        x_dA = (OML.xy[0][iPoint]+OML.xy[0][iPoint-1])/2
        y_dA = (OML.xy[1][iPoint]+OML.xy[1][iPoint-1])/2/2
        IYY_OML_arr[iPoint,0] = abs(x_dA**2*dA[iPoint,0])
        IXY_OML_arr[iPoint,0] = x_dA*y_dA*dA[iPoint,0]
        x_bar_num_arr[iPoint,0] = dA[iPoint,0]*x_dA
    
    A = np.abs(np.sum(dA))
    IXX_OML = np.sum(IXX_OML_arr)
    IYY_OML = np.sum(IYY_OML_arr)
    IXY_OML = np.sum(IXY_OML_arr)
    x_bar = np.sum(x_bar_num_arr)/A
    y_bar = np.sum(y_bar_num_arr)/A
    
    return A, IXX_OML, IYY_OML, IXY_OML, x_bar, y_bar


def constantThickApproxAirfoilProps(profile, scalingFactor, thickness):
    """
    Computes seciton properties of a basic profile of uniform thickness, 
    such as an airfoil. Note that the moment of area and centroids are not 
    properly shifted, an average betweeen IML and OML is assumed for brevity. 

    Parameters
    ----------
    profile : 2D np.array
        Numpy array of profile such as Selig format airfoil data
    scalingFactor : float
        Value for scaling 
    thickness : float
        Value defining wall thickness

    Returns
    -------
    Section properties: 
        Dictionary of section properties: Ixx, Iyy, Ixy, A, x_bar, y_bar

    """
    
    #Use Shapely library to manipulate airfoil profiles and build section
    OML = geometry.polygon.LineString(profile*scalingFactor)
    OML_list = [[x,y] for x,y in zip(OML.xy[0],OML.xy[1])]
    IML = OML.parallel_offset(-thickness)    
    IML_list = [[x,y] for x,y in zip(IML.xy[0],IML.xy[1])]
    
    #OML properties
    OML_A = drangSectionProps.area(OML_list)
    OML_x_bar, OML_y_bar = drangSectionProps.centroid(OML_list)
    OML_Ixx_c, OML_Iyy_c, OML_Ixy_c = drangSectionProps.inertia(OML_list)
    
    #IML properties
    IML_A = drangSectionProps.area(IML_list)
    IML_x_bar, IML_y_bar = drangSectionProps.centroid(IML_list)
    IML_Ixx_c, IML_Iyy_c, IML_Ixy_c = drangSectionProps.inertia(IML_list)
    
    #Net sectin properties - NOTE that the moments of inerta are simply subtracted and the net centroid 
    #is approximated as the average for previty
    sectionProps = {}
    sectionProps['A'] = OML_A - IML_A
    sectionProps['x_bar'] = (OML_x_bar + IML_x_bar)/2
    sectionProps['y_bar'] = (OML_y_bar + IML_y_bar)/2
    sectionProps['Ixx_c'] = OML_Ixx_c-IML_Ixx_c
    sectionProps['Iyy_c'] = OML_Iyy_c-IML_Iyy_c
    sectionProps['Ixy_c'] = OML_Ixy_c-IML_Ixy_c
    
    fig, axs = plt.subplots(1, 1)
    # axs.plot(profile[:,0],profile[:,1],label='Original')
    axs.plot(OML.xy[0], OML.xy[1],label='OML')
    axs.plot(IML.xy[0], IML.xy[1],label='IML')
    axs.set_aspect('equal', 'box')
    axs.set_title('Scaled Airfoil OML/IML')
    axs.set_xlabel('x [m]')
    axs.set_ylabel('y [m]')
    axs.grid()
    axs.legend()
    
    return sectionProps
    