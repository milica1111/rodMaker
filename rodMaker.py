""" 
################################################################
            BlockMeshDict Creation Script for Fuel Rod
===============================================================
################################################################
This script generates the 'blockMeshDict' file for meshing a 1D,
2D-smeared, 2D-discrete, and 3D nuclear fuel rod composed of 
discrete pellets. It automates the  process by reading the necessary
input parameters from a file named 'rodDict'.The 'rodDict' file must
be located in the same directory as this script.
"""

import math
from collections import defaultdict
import os
import re
# importing the module 
import ast
import random

###################################################################################################################################################
######################----------------------------- UNIVERSAL FUNCTIONS FOR ALL GEOMETRIES ------------------------------######################
###################################################################################################################################################

def appendBlock(list_blocks, vertices, mesh, name):
    block = {
        "name" : name,
        "vertices" : vertices,
        "mesh" : mesh
    }
    list_blocks.append(block)

def addToPatchDict(patchDict, name, type, neighbour, owner, face):
    if name not in patchDict:
        patchDict[name] = {
            "type": type,
            "neighbour": neighbour,
            "owner": owner,
            "faces": []
        }
    # this is for the cases when we still do not add any face but
    # add values to other keys
    if face!=[]:
        patchDict[name]["faces"].append(face)


###################################################################################################################################################
######################----------------------------- FUNCTIONS FOR 1D AND 2D-SMEARED GEOMETRIES ------------------------------######################
###################################################################################################################################################

def addWedgeVertices(list_vertices, block, wedgeAngle, offset, clad_flag):
        
    rInner=block["rInner"]
    rOuter=block["rOuter"]
    h=block["height"]
    theta=wedgeAngle

    # due to having just one slice of the rod, the implementation of
    # the arc is useless and the straight faces are used. Vertices are
    # slightly shifted to maintain the total volume of the rod slice

    correction=math.sqrt(theta/math.sin(theta))

    xInner=rInner*math.cos(theta/2)*correction
    yInner=rInner*math.sin(theta/2)*correction
    xOuter=rOuter*math.cos(theta/2)*correction
    yOuter=rOuter*math.sin(theta/2)*correction

    list_vertices.append([xInner,-yInner,offset])
    list_vertices.append([xInner,+yInner,offset])
    list_vertices.append([xOuter,-yOuter,offset])
    list_vertices.append([xOuter,+yOuter,offset])
    list_vertices.append([xInner,-yInner,offset + h])
    list_vertices.append([xInner,+yInner,offset + h])
    list_vertices.append([xOuter,-yOuter,offset + h])
    list_vertices.append([xOuter,+yOuter,offset + h])

    if clad_flag==1:
        type=block["type"]
        if type=="cap":
            list_vertices.append([0.0, 0.0 ,offset])
            #list_vertices.append([0.0, 0.0 ,offset])
            list_vertices.append([0.0, 0.0 ,offset+h])
            #list_vertices.append([0.0, 0.0 ,offset+h])

def addWedgeBlocks(list_blocks, block, i_vertex, clad_flag):
    if clad_flag==0:
        shift=block["nVertices"]/2
    else:
        shift=4

    nR=block["nR"]
    nZ=block["nZ"]
    name=block["name"]
    
    vertices=[x + i_vertex for x in [0, 2, 3, 1]]
    vertices.extend([x + shift for x in vertices])
    appendBlock(list_blocks, vertices, [nR, 1, nZ], name)

    if clad_flag==1:
        type=block["type"]
        if type=="cap":
            nR=block["nRInner"]
            vertices=[x + i_vertex for x in [8, 0, 1, 8, 9, 4, 5, 9]]
            appendBlock(list_blocks, vertices, [nR, 1, nZ], name)

def addWedgePatches(patchDict, mergePatchDict, block, nBlocks, bottomCap, topCap, i_vertex, i_global, geometry, flag):
    shift=4

    if flag==0: # fuel

        rInner=block["rInner"]

        ####################################### Fuel bottom patch #############################################
        base=[x + i_vertex for x in [0, 1, 3, 2]]
        if i_global==1:
                addToPatchDict(patchDict, "fuelBottom", "patch", "none", "false", base)
                if geometry=='1D':
                    patchDict["fuelBottom"]["type"]='empty'
                elif bottomCap:
                    patchDict["fuelBottom"]["type"]='regionCoupledOFFBEAT'
                    patchDict["fuelBottom"]["neighbour"]='bottomCapInner'
                    patchDict["fuelBottom"]["owner"]="true"
        else:
            addToPatchDict(patchDict,  "fuelBottom_" + str(i_global), "patch", "none", "false", base)
        
        ####################################### Fuel top patch #################################################
        base=[x + shift for x in base[::-1]]
        if i_global==nBlocks:
            addToPatchDict(patchDict, "fuelTop", "patch", "none", "true", base)
            if geometry=='1D':
                patchDict["fuelTop"]["type"]='empty'
            elif topCap:
                    patchDict["fuelTop"]["type"]='regionCoupledOFFBEAT'
                    patchDict["fuelTop"]["neighbour"]='topCapInner'
                    patchDict["fuelTop"]["owner"]="false"
        else:
            addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "true", base)
            mergePatchDict["fuelTop_" + str(i_global)] = 'fuelBottom_' + str(i_global+1)
        
        ################################# Fuel front and back patches ##########################################
        base=[x + i_vertex for x in [0, 2, 6, 4]]
        addToPatchDict(patchDict, "fuelFront", "wedge", "none", "false", base)
        base=[x + i_vertex for x in [1, 5, 7, 3]]
        addToPatchDict(patchDict, "fuelBack", "wedge", "none", "false", base)

        ################################ Fuel inner and outer patches ##########################################
        if rInner>0:
            base=[x + i_vertex for x in [0, 4, 5, 1]]
            addToPatchDict(patchDict, "fuelInner", "patch", "none", "false", base)
        
        base=[x + i_vertex for x in [2, 3, 7, 6]]
        addToPatchDict(patchDict, "fuelOuter", "regionCoupledOFFBEAT", "cladInner", "true", base)

    else: # cladding
        type=block["type"]
        ####################################### Cladding bottom patch ##########################################
        base=[x + i_vertex for x in [0, 1, 3, 2]]
        if i_global==1:
            addToPatchDict(patchDict, "cladBottom", "patch", "none", "false", base)
            if bottomCap:
                base=[x + i_vertex for x in [8,1, 0, 8]]
                addToPatchDict(patchDict, "cladBottom", "patch", "none", "false", base)
                base=[x + i_vertex for x in [9, 4, 5, 9]]
                addToPatchDict(patchDict, "bottomCapInner", "regionCoupledOFFBEAT", "fuelBottom", "false", base)
        else:
            addToPatchDict(patchDict, "cladBottom_" + str(i_global), "patch", "none", "false", base)

        ####################################### Cladding top patch #############################################
        base=[x + i_vertex+shift for x in [0, 2, 3, 1]]
        if i_global==nBlocks:
            addToPatchDict(patchDict, "cladTop", "patch", "none", "false", base)
            if topCap:
                base=[x + i_vertex for x in [9, 4, 5, 9]]
                addToPatchDict(patchDict, "cladTop", "patch", "none", "false", base)
                base=[x + i_vertex for x in [8, 1, 0, 8]]
                addToPatchDict(patchDict, "topCapInner", "regionCoupledOFFBEAT", "fuelTop", "true", base)
        else:
            addToPatchDict(patchDict, "cladTop_" + str(i_global), "patch", "none", "false", base)
            mergePatchDict["cladTop_" + str(i_global)] = 'cladBottom_' + str(i_global+1)

        ####################################### Cladding side patches ###########################################
        base=[x + i_vertex for x in [0, 2, 6, 4]]
        addToPatchDict(patchDict, "cladFront", "wedge", "none", "false", base)
        base=[x + i_vertex for x in [1, 5, 7, 3]]
        addToPatchDict(patchDict, "cladBack", "wedge", "none", "false", base)
        if type=="cap":
                addToPatchDict(patchDict, "cladFront", "wedge", "none", "false", [x + i_vertex for x in [8, 0, 4, 9]])
                addToPatchDict(patchDict, "cladBack", "wedge", "none", "false", [x + i_vertex for x in [8, 9, 5, 1]])

        ################################ Cladding inner and outer patches ########################################
        if type=="normal":
            base=[x + i_vertex for x in [0, 4, 5, 1]]
            addToPatchDict(patchDict, "cladInner", "regionCoupledOFFBEAT", "fuelOuter", "false", base)
        
        base=[x + i_vertex for x in [2, 3, 7, 6]]
        addToPatchDict(patchDict, "cladOuter", "patch", "none", "false", base)
                
###################################################################################################################################################
######################----------------------------- FUNCTIONS FOR 3D and 2D-discrete geometries -----------------------------######################
###################################################################################################################################################

def addSpheres(list_spheres,pellet, offset, geometry, shiftX=0, shiftY=0):
    # Calculate the cennter of the sphere
    if geometry=='3D':
        if pellet["type"]=='dished' or pellet["type"]=='dishedChamfered':

            R=pellet["rCurvatureDish"]
            r_dish=pellet["rDish"]
            h=pellet["height"]

            z_bottom=offset - math.sqrt(R**2-r_dish**2)
            z_top=offset + h + math.sqrt(R**2-r_dish**2)

            bottom_sphere = {

                "z":                z_bottom,
                "radius":           R,
                "x":                shiftX,
                "y":                shiftY

            }
            list_spheres.append(bottom_sphere)

            top_sphere = {
                "z":                z_top,
                "radius":           R,
                "x":                shiftX,
                "y":                shiftY
            }
            
            list_spheres.append(top_sphere)


def append4SymVertices(list_vertices, xy, height, shiftX=0, shiftY=0):
    list_vertices.append([-xy+shiftX,-xy+shiftY,height])
    list_vertices.append([+xy+shiftX,-xy+shiftY,height])
    list_vertices.append([+xy+shiftX,+xy+shiftY,height])
    list_vertices.append([-xy+shiftX,+xy+shiftY,height])

def addPelletVertices(list_vertices, pellet, offset, geometry, shiftX=0, shiftY=0):

    type=pellet["type"]
    rInner=pellet["rInner"]
    rOuter=pellet["rOuter"]
    h=pellet["height"]


    #********************************************************************************************
    #**********************************  2D-discrete ********************************************
    #********************************************************************************************

    if geometry=='2D-discrete':

        theta=pellet["wedgeAngle"]
        xInner=rInner*math.cos(theta/2)
        yInner=rInner*math.sin(theta/2)

    ####################################################################
    ######## CASE 1: Dished and chamfered ##############################
    ####################################################################
        
        if type =="dishedChamfered":

            R=pellet["rCurvatureDish"]
            rDish=pellet["rDish"]
            rLand=pellet["rLand"]

            hInner= math.sqrt(R**2 - rInner**2) - math.sqrt(R**2 - rDish**2)

            xDishEnd=rDish*math.cos(theta/2)
            yDishEnd=rDish*math.sin(theta/2)

            xLandEnd=rLand*math.cos(theta/2)
            yLandEnd=rLand*math.sin(theta/2)

            xChamferEnd=rOuter*math.cos(theta/2)
            yChamferEnd=rOuter*math.sin(theta/2)
            hChamfer=pellet["chamferHeight"]

            list_vertices.append([xInner, -yInner, offset + hInner])
            list_vertices.append([xInner, yInner, offset + hInner])
            list_vertices.append([xDishEnd, -yDishEnd, offset])
            list_vertices.append([xDishEnd, yDishEnd, offset])
            list_vertices.append([xLandEnd, -yLandEnd, offset])
            list_vertices.append([xLandEnd, yLandEnd, offset])
            list_vertices.append([xChamferEnd, -yChamferEnd, offset + hChamfer])
            list_vertices.append([xChamferEnd, yChamferEnd, offset + hChamfer])

            list_vertices.append([xInner, -yInner, offset + h - hInner])
            list_vertices.append([xInner, yInner, offset + h - hInner])
            list_vertices.append([xDishEnd, -yDishEnd, offset + h])
            list_vertices.append([xDishEnd, yDishEnd, offset + h])
            list_vertices.append([xLandEnd, -yLandEnd, offset + h])
            list_vertices.append([xLandEnd, yLandEnd, offset + h])
            list_vertices.append([xChamferEnd, -yChamferEnd, offset + h - hChamfer])
            list_vertices.append([xChamferEnd, yChamferEnd, offset + h - hChamfer])

        ###########################################################################
        ################# CASE 2: Just dished #####################################
        ###########################################################################
        elif type=="dished":

            R=pellet["rCurvatureDish"]
            rDish=pellet["rDish"]
        
            hInner= math.sqrt(R**2 - rInner**2) - math.sqrt(R**2 - rDish**2)

            xDishEnd=rDish*math.cos(theta/2)
            yDishEnd=rDish*math.sin(theta/2)

            xEnd=rOuter*math.cos(theta/2)
            yEnd=rOuter*math.sin(theta/2)

            list_vertices.append([xInner, -yInner, offset + hInner])
            list_vertices.append([xInner, yInner, offset + hInner])
            list_vertices.append([xDishEnd, -yDishEnd, offset])
            list_vertices.append([xDishEnd, yDishEnd, offset])
            list_vertices.append([xEnd, -yEnd, offset])
            list_vertices.append([xEnd, yEnd, offset])

            list_vertices.append([xInner, -yInner, offset + h - hInner])
            list_vertices.append([xInner, yInner, offset + h - hInner])
            list_vertices.append([xDishEnd, -yDishEnd, offset + h])
            list_vertices.append([xDishEnd, yDishEnd, offset + h])
            list_vertices.append([xEnd, -yEnd, offset + h])
            list_vertices.append([xEnd, yEnd, offset + h])

        ###########################################################################
        ################# CASE 3: Just chamfered ##################################
        ###########################################################################
        elif type=="chamfered":

            rChamferStart=pellet["rLand"]
            hChamfer=pellet["chamferHeight"]

            xChamferStart=rChamferStart*math.cos(theta/2)
            yChamferStart=rChamferStart*math.sin(theta/2)

            xChamferEnd=rOuter*math.cos(theta/2)
            yChamferEnd=rOuter*math.sin(theta/2)

            list_vertices.append([xInner, -yInner, offset])
            list_vertices.append([xInner, yInner, offset])
            list_vertices.append([xChamferStart, -yChamferStart, offset])
            list_vertices.append([xChamferStart, yChamferStart, offset])
            list_vertices.append([xChamferEnd, -yChamferEnd, offset + hChamfer])
            list_vertices.append([xChamferEnd, yChamferEnd, offset + hChamfer])

            list_vertices.append([xInner, -yInner, offset + h])
            list_vertices.append([xInner, yInner, offset + h])
            list_vertices.append([xChamferStart, -yChamferStart, offset + h])
            list_vertices.append([xChamferStart, yChamferStart, offset + h])
            list_vertices.append([xChamferEnd, -yChamferEnd, offset + h - hChamfer])
            list_vertices.append([xChamferEnd, yChamferEnd, offset + h - hChamfer])

        ###########################################################################
        ################# CASE 4: Flat ############################################
        ########################################################################### 
        else:
            xEnd=rOuter*math.cos(theta/2)
            yEnd=rOuter*math.sin(theta/2)
            list_vertices.append([xInner, -yInner, offset])
            list_vertices.append([xInner, yInner, offset])
            list_vertices.append([xEnd, -yEnd, offset])
            list_vertices.append([xEnd, yEnd, offset])

            list_vertices.append([xInner, -yInner, offset + h])
            list_vertices.append([xInner, yInner, offset + h]) 
            list_vertices.append([xEnd, -yEnd, offset + h])
            list_vertices.append([xEnd, yEnd, offset + h])

    #********************************************************************************************
    #************************************* 3D geometry ******************************************
    #********************************************************************************************
            
    if geometry=='3D':

    ####################################################################
    ######## CASE 1: Dished and chamfered ##############################
    ####################################################################

        if type =="dishedChamfered":

            R=pellet["rCurvatureDish"]
            rDish=pellet["rDish"]
           
            rLand=pellet["rLand"]

            xyDishEnd=rDish/math.sqrt(2)
            xyLandEnd=rLand/math.sqrt(2)
            # hDishEnd=0, hLandEnd=0

            xyChamferEnd=rOuter/math.sqrt(2)
            hChamfer=pellet["chamferHeight"]

            if rInner==0:

                f=pellet['squareFraction']
                rCorner=f*rDish

                xyFirst=rCorner/math.sqrt(2)
                hFirst=math.sqrt(R**2 - rCorner**2) - math.sqrt(R**2 - rDish**2)

            else: #there is central hole

                xyFirst=rInner/math.sqrt(2)
                hFirst= math.sqrt(R**2 - rInner**2) - math.sqrt(R**2 - rDish**2)

            ################# Appending vertices to the list ######################
                
            ################# vertices at the pellet bottom #######################
            append4SymVertices(list_vertices, xyFirst, offset + hFirst, shiftX, shiftY)
            append4SymVertices(list_vertices, xyDishEnd, offset, shiftX, shiftY)
            append4SymVertices(list_vertices, xyLandEnd, offset, shiftX, shiftY)
            append4SymVertices(list_vertices, xyChamferEnd, offset + hChamfer, shiftX, shiftY)

            ################# vertices at the pellet top ##########################
            append4SymVertices(list_vertices, xyFirst, offset + h - hFirst, shiftX, shiftY)
            append4SymVertices(list_vertices, xyDishEnd, offset + h, shiftX, shiftY)
            append4SymVertices(list_vertices, xyLandEnd, offset + h, shiftX, shiftY)
            append4SymVertices(list_vertices, xyChamferEnd, offset + h - hChamfer, shiftX, shiftY)

        ###########################################################################
        ################# CASE 2: Just dished #####################################
        ###########################################################################
            
        elif type=="dished":

            R=pellet["rCurvatureDish"]
            rDish=pellet["rDish"]

            xyDishEnd=rDish/math.sqrt(2)
            # hDishEnd=0

            xyEnd=rOuter/math.sqrt(2)
            # hEnd=0

            if rInner==0:

                f=pellet["squareFraction"]
                rCorner=f*rDish

                xyFirst=rCorner/math.sqrt(2)
                hFirst=math.sqrt(R**2 - rCorner**2) - math.sqrt(R**2 - rDish**2)

            else: #there is central hole

                xyFirst=rInner/math.sqrt(2)
                # * May be the error message added if the inner radius is larger than the dish radius itself
                hFirst= math.sqrt(R**2 - rInner**2) - math.sqrt(R**2 - rDish**2)
            
            ################# Appending vertices to the list ######################
                
            ################# vertices at the pellet bottom #######################
            append4SymVertices(list_vertices, xyFirst, offset + hFirst, shiftX, shiftY)
            append4SymVertices(list_vertices, xyDishEnd, offset, shiftX, shiftY)
            append4SymVertices(list_vertices, xyEnd, offset, shiftX, shiftY)

            ################# vertices at the pellet top ##########################
            append4SymVertices(list_vertices, xyFirst, offset + h - hFirst, shiftX, shiftY)
            append4SymVertices(list_vertices, xyDishEnd, offset + h, shiftX, shiftY)
            append4SymVertices(list_vertices, xyEnd, offset + h, shiftX, shiftY)


        ###########################################################################
        ################# CASE 3: Just chamfered ##################################
        ###########################################################################
            
        elif type=="chamfered":

            rChamferStart=pellet["rLand"]

            xyChamferStart=rChamferStart/math.sqrt(2)
            hChamfer=pellet["chamferHeight"]

            xyEnd=rOuter/math.sqrt(2)
            # hEnd=0

            if rInner==0:

                f=pellet["squareFraction"]
                rCorner=f*rOuter

                xyFirst=rCorner/math.sqrt(2)
                # hFirst=0.0

            else: #there is central hole

                xyFirst=rInner/math.sqrt(2)
                # hFirst=0.0

            ################# Appending vertices to the list ######################
                
            ################# vertices at the pellet bottom #######################
            append4SymVertices(list_vertices, xyFirst, offset, shiftX, shiftY)
            append4SymVertices(list_vertices, xyChamferStart, offset, shiftX, shiftY)
            append4SymVertices(list_vertices, xyEnd, offset + hChamfer, shiftX, shiftY)

            ################# vertices at the pellet top ##########################
            append4SymVertices(list_vertices, xyFirst, offset + h, shiftX, shiftY)
            append4SymVertices(list_vertices, xyChamferStart, offset + h, shiftX, shiftY)
            append4SymVertices(list_vertices, xyEnd, offset + h - hChamfer, shiftX, shiftY)

        ###########################################################################
        ################# CASE 4: Flat ############################################
        ###########################################################################
        else:
            
            xyEnd=rOuter/math.sqrt(2)
            # hEnd=0

            if rInner==0:

                f=pellet["squareFraction"]
                rCorner=f*rOuter

                xyFirst=rCorner/math.sqrt(2)
                # hFirst=0.0

            else: #there is central hole

                xyFirst=rInner/math.sqrt(2)
                # hFirst=0.0

            ################# Appending vertices to the list ######################
                
            ################# vertices at the pellet bottom #######################
            append4SymVertices(list_vertices, xyFirst, offset, shiftX, shiftY)
            append4SymVertices(list_vertices, xyEnd, offset, shiftX, shiftY)

            ################# vertices at the pellet top ##########################
            append4SymVertices(list_vertices, xyFirst, offset + h, shiftX, shiftY)
            append4SymVertices(list_vertices, xyEnd, offset + h, shiftX, shiftY)

def addCladVertices(list_vertices, clad_block, offset, geometry):
    type=clad_block["type"]
    rInner=clad_block["rInner"]
    rOuter=clad_block["rOuter"]
    h=clad_block["height"]

    if geometry=="3D":

        xyInner=rInner/math.sqrt(2)
        xyOuter=rOuter/math.sqrt(2)

        if type=="cap":
            f=clad_block["squareFraction"]
            r=f*rInner
            xySquare=r/math.sqrt(2)
            append4SymVertices(list_vertices, xySquare, offset)

        append4SymVertices(list_vertices, xyInner, offset)
        append4SymVertices(list_vertices, xyOuter, offset)

        if type=="cap":
            append4SymVertices(list_vertices, xySquare, offset+h)

        append4SymVertices(list_vertices, xyInner, offset+h)
        append4SymVertices(list_vertices, xyOuter, offset+h)
    
    else:
        addWedgeVertices(list_vertices, clad_block, clad_block["wedgeAngle"], offset, 1)



def append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, meshXY, meshRadial, meshZ, name):

    vertices=[x + i_vertex for x in baseFace]
    ending=vertices[:2][::-1] #in the last of 4 blocks they will be set to vertices
    vertices.extend([x + shift for x in vertices])
    mesh= [meshRadial, meshXY, meshZ]
    appendBlock(list_blocks, vertices, mesh, name)

    vertices = [x + 1 for x in vertices]
    appendBlock(list_blocks, vertices, mesh, name)

    vertices = [x + 1 for x in vertices]
    appendBlock(list_blocks, vertices, mesh, name)

    vertices = [x + 1 for x in vertices]
    vertices[2:4]=ending
    vertices[-2:]=[x + shift for x in ending]
    appendBlock(list_blocks, vertices, mesh, name)


def addFuelBlocks(list_blocks, pellet, i_vertex, geometry):

    type=pellet["type"]
    rInner=pellet["rInner"]
    shift=pellet["nVertices"]/2
    nCellsZ=pellet["nCellsZPellet"]
    nCellsRTotal=pellet["nCellsRPellet"]
    name=pellet["blockName"]

    if geometry=="2D-discrete":
        ####################################################################
        ######## CASE 1: Dished and chamfered ##############################
        ####################################################################
        if type=="dishedChamfered":

            nCellsRDish=pellet["nCellsRDish"]
            nCellsRChamfer=pellet["nCellsRChamfer"]
            nCellsRLand=nCellsRTotal-nCellsRDish-nCellsRChamfer

            vertices=[x + i_vertex for x in [0, 2, 3, 1]]
            vertices.extend([x + shift for x in vertices])
            mesh= [nCellsRDish, 1 , nCellsZ]
            appendBlock(list_blocks, vertices, mesh, name)

            vertices=[x + 2 for x in vertices]
            mesh= [nCellsRLand, 1 , nCellsZ]
            appendBlock(list_blocks, vertices, mesh, name)

            vertices=[x + 2 for x in vertices]
            mesh= [nCellsRChamfer, 1 , nCellsZ]
            appendBlock(list_blocks, vertices, mesh, name)

        ###########################################################################
        ################# CASE 2: Just dished #####################################
        ###########################################################################
        elif type=="dished":

            nCellsRDish=pellet["nCellsRDish"]
            nCellsRLand=nCellsRTotal-nCellsRDish

            vertices=[x + i_vertex for x in [0, 2, 3, 1]]
            vertices.extend([x + shift for x in vertices])
            mesh= [nCellsRDish, 1 , nCellsZ]
            appendBlock(list_blocks, vertices, mesh, name)

            vertices=[x + 2 for x in vertices]
            mesh= [nCellsRLand, 1 , nCellsZ]
            appendBlock(list_blocks, vertices, mesh, name)

        ###########################################################################
        ################# CASE 3: Just chamfered ##################################
        ###########################################################################
        elif type=="chamfered":

            nCellsRChamfer=pellet["nCellsRChamfer"]
            nCellsRLand=nCellsRTotal-nCellsRChamfer

            vertices=[x + i_vertex for x in [0, 2, 3, 1]]
            vertices.extend([x + shift for x in vertices])
            mesh= [nCellsRLand, 1 , nCellsZ]
            appendBlock(list_blocks, vertices, mesh, name)

            vertices=[x + 2 for x in vertices]
            mesh= [nCellsRChamfer, 1 , nCellsZ]
            appendBlock(list_blocks, vertices, mesh, name)

        ###########################################################################
        ################# CASE 4: Flat ############################################
        ###########################################################################
        else:
            vertices=[x + i_vertex for x in [0, 2, 3, 1]]
            vertices.extend([x + shift for x in vertices])
            mesh= [nCellsRTotal, 1 , nCellsZ]
            appendBlock(list_blocks, vertices, mesh, name)



    elif geometry=="3D":

        nCellsAzimuthal=math.ceil(pellet["nCellsAzimuthal"]/4)
        ####################################################################
        ######## CASE 1: Dished and chamfered ##############################
        ####################################################################
        if type=="dishedChamfered":
                
                nCellsRDish=pellet["nCellsRDish"]
                nCellsRChamfer=pellet["nCellsRChamfer"]
                nCellsRLand=nCellsRTotal-nCellsRDish-nCellsRChamfer
                

                if rInner==0:

                    # The block corresponding to the central square
                    nCellsXYSquare=nCellsAzimuthal
                    vertices=[x + i_vertex for x in [0, 1, 2, 3]]
                    vertices.extend([x + shift for x in vertices])
                    mesh= [nCellsXYSquare, nCellsXYSquare, nCellsZ]
                    appendBlock(list_blocks, vertices, mesh, name)

                    # Mesh od the rest of the dish
                    nCellsRadial=math.ceil(nCellsRDish-nCellsXYSquare)
                    nCellsXY=nCellsXYSquare

                else: # when pellet has the central hole
                    nCellsXY=nCellsAzimuthal
                    nCellsRadial=nCellsRDish

                baseFace=[0, 4, 5, 1]
                append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, nCellsXY, nCellsRadial, nCellsZ, name)
                # Adding blocks for land
                baseFace=[x+4 for x in baseFace]
                append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, nCellsXY, nCellsRLand, nCellsZ, name)
                # Adding blocks for chamferred part of the pellet
                baseFace=[x+4 for x in baseFace]
                append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, nCellsXY, nCellsRChamfer, nCellsZ, name)
        

        ###########################################################################
        ################# CASE 2: Just dished #####################################
        ###########################################################################

        elif type=="dished":
            nCellsRDish=pellet["nCellsRDish"]
            nCellsRLand=nCellsRTotal-nCellsRDish

            if rInner==0:
                nCellsXYSquare=nCellsAzimuthal
                vertices=[x + i_vertex for x in [0, 1, 2, 3]]
                vertices.extend([x + shift for x in vertices])
                mesh= [nCellsXYSquare, nCellsXYSquare, nCellsZ]
                appendBlock(list_blocks, vertices, mesh, name)

                # Mesh od the rest of the dish
                nCellsRadial=nCellsRDish-nCellsXYSquare
                nCellsXY=nCellsXYSquare
            else:
                nCellsXY=nCellsAzimuthal
                nCellsRadial=nCellsRDish
            
            baseFace=[0, 4, 5, 1]
            append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, nCellsXY, nCellsRadial, nCellsZ, name)
            # Adding blocks for land, or the rest of the pellet
            baseFace=[x+4 for x in baseFace]
            append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, nCellsXY, nCellsRLand, nCellsZ, name)


        ###########################################################################
        ################# CASE 3: Just chamfered ##################################
        ###########################################################################
        elif type=="chamfered":

            nCellsRChamfer=pellet["nCellsRChamfer"]
            nCellsRLand=nCellsRTotal-nCellsRChamfer

            if rInner==0:
                nCellsXYSquare=nCellsAzimuthal
                vertices=[x + i_vertex for x in [0, 1, 2, 3]]
                vertices.extend([x + shift for x in vertices])
                mesh= [nCellsXYSquare, nCellsXYSquare, nCellsZ]
                appendBlock(list_blocks, vertices, mesh, name)

                # Mesh od the rest of the pellet up to chamferred part
                nCellsRadial=nCellsRLand-nCellsXYSquare
                nCellsXY=nCellsXYSquare
            else:
                nCellsXY=nCellsAzimuthal
                nCellsRadial=nCellsRLand
            
            baseFace=[0, 4, 5, 1]
            append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, nCellsXY, nCellsRadial, nCellsZ, name)
            # Adding blocks for land, or the rest of the pellet
            baseFace=[x+4 for x in baseFace]
            append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, nCellsXY, nCellsRChamfer, nCellsZ, name)
        
        ###########################################################################
        ################# CASE 4: Flat ############################################
        ###########################################################################
            
        else:

            if rInner==0:
                nCellsXYSquare=nCellsAzimuthal
                vertices=[x + i_vertex for x in [0, 1, 2, 3]]
                vertices.extend([x + shift for x in vertices])
                mesh= [nCellsXYSquare, nCellsXYSquare, nCellsZ]
                appendBlock(list_blocks, vertices, mesh, name)

                # Mesh od the rest of the pellet up to chamferred part
                nCellsRadial=nCellsRTotal-nCellsXYSquare
                nCellsXY=nCellsXYSquare
            else:
                nCellsXY=nCellsAzimuthal
                nCellsRadial=nCellsRTotal

            baseFace=[0, 4, 5, 1]
            append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, nCellsXY, nCellsRadial, nCellsZ, name)

def addCladBlocks(list_blocks, clad_block, i_vertex, geometry):

    type=clad_block["type"]
    nCellsZ=clad_block["nCellsZ"]
    nROuter=clad_block["nCellsR"]
    name=clad_block["blockName"]
    
    if geometry=="3D":
        shift=clad_block["nVertices"]/2
        nAzimuthalOuter=clad_block["nCellsAzimuthal"]/4

        if type=="cap":

            # We need to add the square in the center
            nRInner=clad_block["nCellsRInner"]

            vertices=[x + i_vertex for x in [0, 1, 2, 3]]
            vertices.extend([x + shift for x in vertices])
            mesh= [nAzimuthalOuter, nAzimuthalOuter, nCellsZ]
            appendBlock(list_blocks, vertices, mesh, name)

            baseFace=[0,4,5,1]
            nCellsRRest=nRInner-nAzimuthalOuter
            append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, nAzimuthalOuter, nCellsRRest, nCellsZ, name)

            baseFace=[4,8,9,5]
        else:
            baseFace=[0,4,5,1]

        append4AzimuthallySymmBlocks(list_blocks, i_vertex, baseFace, shift, nAzimuthalOuter, nROuter, nCellsZ, name)
        
    else:
            shift=4
            vertices=[x + i_vertex for x in [0, 2, 3, 1]]
            vertices.extend([x + shift for x in vertices])
            mesh= [nROuter, 1 , nCellsZ]
            appendBlock(list_blocks, vertices, mesh, name)

            if type=="cap":
                nRInner=clad_block["nCellsRInner"]
                vertices=[x + i_vertex for x in [8, 0, 1, 8, 9, 4, 5, 9]]
                mesh= [nRInner, 1, nCellsZ]
                appendBlock(list_blocks, vertices, mesh, name)


def appendEdge(list_edges, vertices, midpoint):

    edge = {
        "vertices": vertices,
        "midpoint": midpoint
    }
    list_edges.append(edge)

def append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, pelletHeight, xy, h, offset, shiftX=0, shiftY=0):
    baseEdge=[x + i_vertex for x in baseEdge]

    for i in range(2):
        appendEdge(list_edges, baseEdge, [shiftX, shiftY-xy, offset + h])
        appendEdge(list_edges, [x + 1 for x in baseEdge], [xy+shiftX, shiftY, offset + h])
        appendEdge(list_edges, [x + 2 for x in baseEdge], [shiftX, xy+shiftY, offset + h])
        last=[x + 3 for x in baseEdge]
        last[1]=baseEdge[0]
        appendEdge(list_edges, last, [shiftX-xy, shiftY, offset + h])
        baseEdge=[x + shift for x in baseEdge]
        h=pelletHeight-h

def append8AzimuthallyShiftedSymmEdges(list_edges, i_vertex, baseEdge, shift, pelletHeight, xy, h, offset, shiftX=0,shiftY=0):
    baseEdge=[x + i_vertex for x in baseEdge]
    for i in range(2):
        appendEdge(list_edges, baseEdge, [shiftX-xy, shiftY-xy, offset + h])
        appendEdge(list_edges, [x + 1 for x in baseEdge], [shiftX+xy, shiftY-xy, offset + h])
        appendEdge(list_edges, [x + 2 for x in baseEdge], [shiftX+xy, shiftY+xy, offset + h])
        appendEdge(list_edges, [x + 3 for x in baseEdge], [shiftX-xy, shiftY+xy, offset + h])
        baseEdge=[x + shift for x in baseEdge]
        h=pelletHeight-h

def addFuelEdges(list_edges, pellet, i_vertex, offset, geometry, shiftX,shiftY):

    type=pellet["type"]
    rInner=pellet["rInner"]
    rOuter=pellet["rOuter"]
    shift=pellet["nVertices"]/2
    height=pellet["height"]

    if geometry=="2D-discrete":

        theta=pellet["wedgeAngle"]/2
        if type=="dishedChamfered" or type=="dished":

            R=pellet["rCurvatureDish"]
            rDish=pellet["rDish"]

            hArcPoint=math.sqrt(R**2-(rDish/2)**2)-math.sqrt(R**2-rDish**2)
            xArcPoint=rDish/2*math.cos(theta/2)
            yArcPoint=rDish/2*math.sin(theta/2)

            vertices1=[x + i_vertex for x in [0, 2]]
            vertices2=[x + 1 for x in vertices1]

            appendEdge(list_edges, vertices1 , [xArcPoint, -yArcPoint, offset + hArcPoint])
            appendEdge(list_edges, vertices2, [xArcPoint, yArcPoint, offset + hArcPoint])
            appendEdge(list_edges, [x + shift for x in vertices1] , [xArcPoint, -yArcPoint, offset+ height - hArcPoint])
            appendEdge(list_edges, [x + shift for x in vertices2] , [xArcPoint, yArcPoint, offset + height - hArcPoint])
    else:

        ####################################################################
        ######## CASE 1: Dished and chamfered ##############################
        ####################################################################
        if type=="dishedChamfered":

            R=pellet["rCurvatureDish"]
            rDish=pellet["rDish"]
            rLand=pellet["rLand"]
            chamferHeight=pellet["chamferHeight"]

            if rInner==0:
                f=pellet["squareFraction"]
                rCorner=f*rDish
                a=rCorner/math.sqrt(2)

                H0=math.sqrt(R**2-rCorner**2)
                # R1^2=R^2-rCorner^2+a^2, R1 is the distance
                # from the center of the sphere to the
                # mid of the edge of the square
                R1=math.sqrt(H0**2+a**2)

                # from triangle similarity:
                H=H0*R/R1

                xy1=a*R/R1
                h1= H - math.sqrt(R**2-rDish**2)

                r2=(1+f)*rDish/2
            else:
                # the middle of the arc:
                xy1=rInner
                h1=math.sqrt(R**2-rInner**2)- math.sqrt(R**2-rDish**2)
                r2=(rDish+rInner)/2

            baseEdge1=[0, 1]
            append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge1, shift, height, xy1, h1, offset,shiftX,shiftY)

            xy2=r2/math.sqrt(2)
            h2=math.sqrt(R**2-r2**2)-math.sqrt(R**2-rDish**2)
            baseEdge2=[0, 4]
            append8AzimuthallyShiftedSymmEdges(list_edges, i_vertex, baseEdge2, shift, height, xy2, h2, offset,shiftX,shiftY)

            # arc of the end of the dished region
            xy3=rDish
            h3=0.0
            baseEdge1=[x + 4 for x in baseEdge1]
            append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge1, shift, height, xy3, h3, offset,shiftX,shiftY)

            # arc of the end of land region
            xy4=rLand
            h4=0.0
            baseEdge1=[x + 4 for x in baseEdge1]
            append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge1, shift, height, xy4, h4, offset,shiftX,shiftY)

            xy5=rOuter
            h5=chamferHeight
            baseEdge1=[x + 4 for x in baseEdge1]
            append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge1, shift, height, xy5, h5, offset,shiftX,shiftY)

        ###########################################################################
        ################# CASE 2: Just dished #####################################
        ###########################################################################
        elif type=="dished":

            R=pellet["rCurvatureDish"]
            rDish=pellet["rDish"]

            if rInner==0:
                f=pellet["squareFraction"]
                rCorner=f*rDish
                a=rCorner/math.sqrt(2)

                H0=math.sqrt(R**2-rCorner**2)
                # R1^2=R^2-rCorner^2+a^2, R1 is the distance
                # from the center of the sphere to the
                # mid of the edge of the square
                R1=math.sqrt(H0**2+a**2)

                # from triangle similarity:
                H=H0*R/R1

                xy1=a*R/R1
                h1= H - math.sqrt(R**2-rDish**2)

                r2=(1+f)*rDish/2
            else:
                # the middle of the arc:
                xy1=rInner
                h1=math.sqrt(R**2-rInner**2)- math.sqrt(R**2-rDish**2)
                r2=(rDish+rInner)/2

            baseEdge=[0, 1]
            append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, height, xy1, h1, offset,shiftX,shiftY)

            xy2=r2/math.sqrt(2)
            h2=math.sqrt(R**2-r2**2)-math.sqrt(R**2-rDish**2)
            baseEdge2=[0, 4]
            append8AzimuthallyShiftedSymmEdges(list_edges, i_vertex, baseEdge2, shift, height, xy2, h2, offset,shiftX,shiftY)

            # arc of the end of the dished region
            xy3=rDish
            h3=0.0
            baseEdge=[x + 4 for x in baseEdge]
            append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, height, xy3, h3, offset,shiftX,shiftY)

            # arc of the end of the pellet
            xy4=rOuter
            h4=0.0
            baseEdge=[x + 4 for x in baseEdge]
            append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, height, xy4, h4, offset,shiftX,shiftY)

        ###########################################################################
        ################# CASE 3: Just chamfered ##################################
        ###########################################################################
        elif type=="chamfered":
            chamferHeight=pellet["chamferHeight"]
            rLand=pellet["rLand"]

            # Since there is no dish, when rInner==0, the straight lines 
            # should be formed between the vertices of square, so no arcs
            # are needed.

            if rInner!=0:
                # the middle of the arc:
                xy1=rInner
                h1=0.0
                baseEdge=[0, 1]
                append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, height, xy1, h1, offset,shiftX,shiftY)
                

            # arc of the end of the land 
            xy2=rLand
            h2=0.0
            baseEdge=[4, 5]
            append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, height, xy2, h2, offset,shiftX,shiftY)

            # arc of the end of the pellet
            xy3=rOuter
            h3=chamferHeight
            baseEdge=[x + 4 for x in baseEdge]
            append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, height, xy3, h3, offset,shiftX,shiftY)
        
        ###########################################################################
        ################# CASE 4: Flat ############################################
        ###########################################################################
            
        elif type=="flat":

            if rInner!=0:
                # the middle of the arc:
                xy1=rInner
                h1=0.0
                baseEdge=[0, 1]
                append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, height, xy1, h1, offset,shiftX,shiftY)

            xy2=rOuter
            h2=0.0
            baseEdge=[4, 5]
            append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, height, xy2, h2, offset,shiftX,shiftY)


def addCladEdges(list_edges, clad_block, i_vertex, offset, geometry):

    if geometry=="3D":
        rInner= clad_block["rInner"]
        rOuter= clad_block["rOuter"]
        height=clad_block["height"]
        type=clad_block["type"]
        shift=clad_block["nVertices"]/2

        if type=="cap":
            baseEdge= [4, 5]
        else:
            baseEdge=[0, 1]

        xy=rInner
        h=0.0
        append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, height, xy, h, offset)
        
        xy=rOuter
        h=0.0

        baseEdge=[x + 4 for x in baseEdge]
        append8AzimuthallySymmEdges(list_edges, i_vertex, baseEdge, shift, height, xy, h, offset)


def appendFaceProjection(list_projection_faces, face, sphere_index):

    projection = {
        "face" : face,
        "sphere" :  f'sphere_{sphere_index}'
    }
    list_projection_faces.append(projection)

def addFaceProjections(list_projection_faces, pellet, i_vertex, i_sphere, geometry):
    if geometry=="3D":
        type=pellet["type"]
        rInner=pellet["rInner"]

        if type=="dishedChamfered" or type=="dished":
            shift=pellet["nVertices"]/2

            # projections on the bottom sphere
            if rInner==0:
                base1=[x + i_vertex for x in [0, 3, 2, 1]]
                appendFaceProjection(list_projection_faces, base1, i_sphere-2)
            base2=[x + i_vertex for x in [4, 0, 1, 5]]
            appendFaceProjection(list_projection_faces, base2, i_sphere-2)
            appendFaceProjection(list_projection_faces, [x + 1 for x in base2], i_sphere-2)
            appendFaceProjection(list_projection_faces, [x + 2 for x in base2], i_sphere-2)
            vector = [x + 3 for x in base2]
            vector[-2:] = base2[:2][::-1]
            appendFaceProjection(list_projection_faces, vector, i_sphere-2)

            # projections on the top sphere
            if rInner==0:
                appendFaceProjection(list_projection_faces, [x + shift for x in base1[::-1]], i_sphere-1)
                
            appendFaceProjection(list_projection_faces, [x + shift for x in base2[::-1]], i_sphere-1)
            appendFaceProjection(list_projection_faces, [x + (shift + 1) for x in base2[::-1]], i_sphere-1)
            appendFaceProjection(list_projection_faces, [x + (shift + 2) for x in base2[::-1]], i_sphere-1)
            appendFaceProjection(list_projection_faces, [x + shift for x in vector[::-1]], i_sphere-1)

            i_sphere+=2

        return i_sphere
    else:
        return 0


def append4SymmetricFacestoPatch(patchDict, name, base, side):  
    patchDict[name]["faces"].append(base)
    patchDict[name]["faces"].append([x + 1 for x in base])
    patchDict[name]["faces"].append([x + 2 for x in base])
    vector = [x + 3 for x in base]
    if side=="bottom":
        vector[-2:] = base[:2][::-1]
    else: # "top"
        vector[:2] = base[-2:][::-1]
    
    patchDict[name]["faces"].append(vector)

    

def addFuelToPatchDict(patchDict, mergePatchDict, pellet, merging, N_pellets, bottomCap, topCap, i_vertex, i_global, geometry):

    type=pellet["type"]
    shift=int(pellet["nVertices"]/2)
    rInner=pellet["rInner"]

#************************************************   2D-discrete  ******************************************************
    if geometry=="2D-discrete":
    ###################################################################################################################
    ####################################  Setting the bottom boundary patches #########################################
    ###################################################################################################################
        base=[x + i_vertex for x in [0, 1, 3, 2]]
        base1=[x + 2 for x in base]
        base2=[x + 2 for x in base1]
        if i_global==1:
            addToPatchDict(patchDict, "fuelBottom", "patch", "none", "false", base)

            if type=="dishedChamfered":
                addToPatchDict(patchDict, "fuelBottom", "patch", "none", "false", base1)
                addToPatchDict(patchDict, "fuelBottom", "patch", "none", "false", base2)
            elif type=="dished" or type=="chamfered":
                addToPatchDict(patchDict, "fuelBottom", "patch", "none", "false", base1)

            if bottomCap:
                patchDict["fuelBottom"]["type"]='regionCoupledOFFBEAT'
                patchDict["fuelBottom"]["neighbour"]='bottomCapInner'
                patchDict["fuelBottom"]["owner"]="true"

        else: # not the bottom pellet
            if type=="dishedChamfered":

                addToPatchDict(patchDict, "dishChamferBottom_" + str(i_global), "patch", "none", "false", base)
                addToPatchDict(patchDict, "dishChamferBottom_" + str(i_global), "patch", "none", "false", base2)
                if not(merging):
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), 'regionCoupledOFFBEAT', "fuelTop_" + str(i_global-1), "false", base1)
                else:
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), "patch", "none", "false", base1)

            elif type=="dished":

                addToPatchDict(patchDict, "dishChamferBottom_" + str(i_global), "patch", "none", "false", base)
                if not(merging):
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), 'regionCoupledOFFBEAT', "fuelTop_" + str(i_global-1), "false", base1)
                else:
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), "patch", "none", "false", base1)

            elif type=="chamfered":

                if not(merging):
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), 'regionCoupledOFFBEAT', "fuelTop_" + str(i_global-1), "false", base)
                else:
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), "patch", "none", "false", base)
                addToPatchDict(patchDict, "dishChamferBottom_" + str(i_global), "patch", "none", "false", base1)

            elif type=="flat":

                if not(merging):
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), 'regionCoupledOFFBEAT', "fuelTop_" + str(i_global-1), "false", base)
                else:
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), "patch", "none", "false", base)


        ###################################################################################################################
        ####################################  Setting the top boundary patches ############################################
        ###################################################################################################################
                
        base=[x + shift for x in base[::-1]]
        base1=[x + shift for x in base1[::-1]]
        base2=[x + shift for x in base2[::-1]]

        if i_global==N_pellets:
            addToPatchDict(patchDict, "fuelTop", "patch", "none", "true", base)
            if type=="dishedChamfered":
                addToPatchDict(patchDict, "fuelTop", "patch", "none", "true", base1)
                addToPatchDict(patchDict, "fuelTop", "patch", "none", "true", base2)
            elif type=="dished" or type=="chamfered":
                addToPatchDict(patchDict, "fuelTop", "patch", "none", "true", base1)            
            if topCap:
                patchDict["fuelTop"]["type"]='regionCoupledOFFBEAT'
                patchDict["fuelTop"]["neighbour"]='topCapInner'
                patchDict["fuelTop"]["owner"]="false"

        else: # not the top pellet

            if type=="dishedChamfered":

                addToPatchDict(patchDict, "dishChamferTop_" + str(i_global), "patch", "none", "false", base)
                addToPatchDict(patchDict, "dishChamferTop_" + str(i_global), "patch", "none", "false", base2)

                if not(merging):
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), 'regionCoupledOFFBEAT', "fuelBottom_" + str(i_global+1), "true", base1)
                else:
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "true", base1)
                    mergePatchDict["fuelTop_" + str(i_global)] = 'fuelBottom_' + str(i_global+1)

            elif type=="dished":

                addToPatchDict(patchDict, "dishChamferTop_" + str(i_global), "patch", "none", "false", base)

                if not(merging):
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), 'regionCoupledOFFBEAT', "fuelBottom_" + str(i_global+1), "true", base1)
                else:
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "true", base1)
                    mergePatchDict["fuelTop_" + str(i_global)] = 'fuelBottom_' + str(i_global+1)

            elif type=="chamfered":

                if not(merging):
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), 'regionCoupledOFFBEAT', "fuelBottom_" + str(i_global-1), "true", base)
                else:
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "true", base)
                    mergePatchDict["fuelTop_" + str(i_global)] = 'fuelBottom_' + str(i_global+1)

                addToPatchDict(patchDict, "dishChamferTop_" + str(i_global), "patch", "none", "false", base1)

            elif type=="flat":

                if not(merging):
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), 'regionCoupledOFFBEAT', "fuelBottom_" + str(i_global+1), "false", base)
                else:
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "false", base)
        
            ###################################################################################################################
            ####################################  Setting the inner boundary patches ##########################################
            ###################################################################################################################
                    
            if rInner>0:
                base=[x + i_vertex for x in [1, 0, shift, shift+1]]
                addToPatchDict(patchDict, "fuelInner", "patch", "none", "false", base)
            
            ###################################################################################################################
            ####################################  Setting the outer boundary patches ##########################################
            ###################################################################################################################

            base=[x + i_vertex for x in [shift-2, shift-1, 2*shift-2, 2*shift-1]]  
            addToPatchDict(patchDict, "fuelOuter", "regionCoupledOFFBEAT", "cladInner", "true", base)
        

    #**************************************************   3D  *********************************************************
    ###################################################################################################################
    ####################################  Setting the bottom boundary patches #########################################
    ###################################################################################################################
    elif geometry=="3D": # 3D geometry
        base=[x + i_vertex for x in [0, 3, 2, 1]]
        base1=[x + i_vertex for x in [4, 0, 1, 5]]
        base2=[x + 4 for x in base1]
        base3=[x + 4 for x in base2]

        if i_global==1:

            addToPatchDict(patchDict, "fuelBottom", "patch", "none", "false", [])
            if rInner==0:
                addToPatchDict(patchDict, "fuelBottom", "patch", "none", "false", base)

            if type=="dishedChamfered":
                append4SymmetricFacestoPatch(patchDict, "fuelBottom", base1, "bottom")
                append4SymmetricFacestoPatch(patchDict, "fuelBottom", base2, "bottom")
                append4SymmetricFacestoPatch(patchDict, "fuelBottom", base3, "bottom")
            elif type=="dished" or type=="chamfered":
                append4SymmetricFacestoPatch(patchDict, "fuelBottom", base1, "bottom")
                append4SymmetricFacestoPatch(patchDict, "fuelBottom", base2, "bottom")
            elif type=="flat":
                append4SymmetricFacestoPatch(patchDict, "fuelBottom", base1, "bottom")
            
            if bottomCap:
                patchDict["fuelBottom"]["type"]='regionCoupledOFFBEAT'
                patchDict["fuelBottom"]["neighbour"]='bottomCapInner'
                patchDict["fuelBottom"]["owner"]="true"
            
        else: #not the bottom pellet

            if type=="dishedChamfered":

                addToPatchDict(patchDict, "dishChamferBottom_" + str(i_global), "patch", "none", "false", [])

                if rInner==0:
                    addToPatchDict(patchDict, "dishChamferBottom_" + str(i_global), "patch", "none", "false", base)

                append4SymmetricFacestoPatch(patchDict, "dishChamferBottom_" + str(i_global), base1, "bottom")
                append4SymmetricFacestoPatch(patchDict, "dishChamferBottom_" + str(i_global), base3, "bottom")

                if not(merging):
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), 'regionCoupledOFFBEAT', "fuelTop_" + str(i_global-1), "false", [])
                else:
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), "patch", "none", "false", [])

                append4SymmetricFacestoPatch(patchDict, "fuelBottom_" + str(i_global), base2, "bottom")

            elif type=="dished":

                addToPatchDict(patchDict, "dishChamferBottom_" + str(i_global), "patch", "none", "false", [])
                if rInner==0:
                    addToPatchDict(patchDict, "dishChamferBottom_" + str(i_global), "patch", "none", "false", base)

                append4SymmetricFacestoPatch(patchDict, "dishChamferBottom_" + str(i_global), base1, "bottom")

                if not(merging):
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), 'regionCoupledOFFBEAT', "fuelTop_" + str(i_global-1), "false", [])
                else:
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), "patch", "none", "false", [])

                append4SymmetricFacestoPatch(patchDict, "fuelBottom_" + str(i_global), base2, "bottom")

            elif type=="chamfered":

                if not(merging):
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), 'regionCoupledOFFBEAT', "fuelTop_" + str(i_global-1), "false", [])
                    if rInner==0:
                        addToPatchDict(patchDict, "fuelBottom_" + str(i_global), 'regionCoupledOFFBEAT', "fuelTop_" + str(i_global-1), "false", base)
                else:
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), "patch", "none", "false", [])
                    if rInner==0:
                        addToPatchDict(patchDict, "fuelBottom_" + str(i_global), "patch", "none", "false", base)

                append4SymmetricFacestoPatch(patchDict, "fuelBottom_" + str(i_global), base1, "bottom")

                addToPatchDict(patchDict, "dishChamferBottom_" + str(i_global), "patch", "none", "false", [])
                append4SymmetricFacestoPatch(patchDict, "dishChamferBottom_" + str(i_global), base2, "bottom")

            elif type=="flat":

                if not(merging):
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), 'regionCoupledOFFBEAT', "fuelTop_" + str(i_global-1), "false", [])
                    if rInner==0:
                        addToPatchDict(patchDict, "fuelBottom_" + str(i_global), 'regionCoupledOFFBEAT', "fuelTop_" + str(i_global-1), "false", base)
                else:
                    addToPatchDict(patchDict, "fuelBottom_" + str(i_global), "patch", "none", "false", [])
                    if rInner==0:
                        addToPatchDict(patchDict, "fuelBottom_" + str(i_global), "patch", "none", "false", base)

                append4SymmetricFacestoPatch(patchDict, "fuelBottom_" + str(i_global), base1, "bottom")

        ###################################################################################################################
        ####################################  Setting the top boundary patches ############################################
        ###################################################################################################################
                
        base=[x + shift for x in base[::-1]]
        base1=[x + shift for x in base1[::-1]]
        base2=[x + shift for x in base2[::-1]]
        base3=[x + shift for x in base3[::-1]]

        if i_global==N_pellets:
                
            addToPatchDict(patchDict, "fuelTop", "patch", "none", "true", [])
            if rInner==0:
                addToPatchDict(patchDict, "fuelTop", "patch", "none", "true", base)

            if type=="dishedChamfered":
                append4SymmetricFacestoPatch(patchDict, "fuelTop", base1, "top")
                append4SymmetricFacestoPatch(patchDict, "fuelTop", base2, "top")
                append4SymmetricFacestoPatch(patchDict, "fuelTop", base3, "top")
            elif type=="dished" or type=="chamfered":
                append4SymmetricFacestoPatch(patchDict, "fuelTop", base1, "top")
                append4SymmetricFacestoPatch(patchDict, "fuelTop", base2, "top")
            elif type=="flat":
                append4SymmetricFacestoPatch(patchDict, "fuelTop", base1, "top")
            if topCap:
                patchDict["fuelTop"]["type"]='regionCoupledOFFBEAT'
                patchDict["fuelTop"]["neighbour"]='topCapInner'
                patchDict["fuelTop"]["owner"]="false"

        else: # not the top pellet

            if type=="dishedChamfered":

                addToPatchDict(patchDict, "dishChamferTop_" + str(i_global), "patch", "none", "false", [])
                if rInner==0:
                    addToPatchDict(patchDict, "dishChamferTop_" + str(i_global), "patch", "none", "false", base)

                append4SymmetricFacestoPatch(patchDict, "dishChamferTop_" + str(i_global), base1, "top")
                append4SymmetricFacestoPatch(patchDict, "dishChamferTop_" + str(i_global), base3, "top")

                if not(merging):
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), 'regionCoupledOFFBEAT', "fuelBottom_" + str(i_global+1), "true", [])
                else:
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "true", [])
                    mergePatchDict["fuelTop_" + str(i_global)] = 'fuelBottom_' + str(i_global+1)

                append4SymmetricFacestoPatch(patchDict, "fuelTop_" + str(i_global), base2, "top")

            elif type=="dished":

                addToPatchDict(patchDict, "dishChamferTop_" + str(i_global), "patch", "none", "false", [])
                if rInner==0:
                    addToPatchDict(patchDict, "dishChamferTop_" + str(i_global), "patch", "none", "false", base)

                append4SymmetricFacestoPatch(patchDict, "dishChamferTop_" + str(i_global), base1, "top")

                if not(merging):
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), 'regionCoupledOFFBEAT', "fuelBottom_" + str(i_global+1), "true", [])
                else:
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "true", [])
                    mergePatchDict["fuelTop_" + str(i_global)] = 'fuelBottom_' + str(i_global+1)

                append4SymmetricFacestoPatch(patchDict, "fuelTop_" + str(i_global), base2, "top")

            elif type=="chamfered":

                if not(merging):
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), 'regionCoupledOFFBEAT', "fuelBottom_" + str(i_global-1), "true", [])
                    if rInner==0:
                        addToPatchDict(patchDict, "fuelTop_" + str(i_global), 'regionCoupledOFFBEAT', "fuelBottom_" + str(i_global-1), "true", base)
                else:
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "true", [])
                    mergePatchDict["fuelTop_" + str(i_global)] = 'fuelBottom_' + str(i_global+1)
                    if rInner==0:
                        addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "true", base)

                append4SymmetricFacestoPatch(patchDict, "fuelTop_" + str(i_global), base1, "top")

                addToPatchDict(patchDict, "dishChamferTop_" + str(i_global), "patch", "none", "false", [])
                append4SymmetricFacestoPatch(patchDict, "dishChamferTop_" + str(i_global), base2, "top")

            elif type=="flat":

                if not(merging):
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), 'regionCoupledOFFBEAT', "fuelBottom_" + str(i_global+1), "false", [])
                    if rInner==0:
                        addToPatchDict(patchDict, "fuelTop_" + str(i_global), 'regionCoupledOFFBEAT', "fuelBottom_" + str(i_global+1), "false", base)
                else:
                    addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "false", [])
                    if rInner==0:
                        addToPatchDict(patchDict, "fuelTop_" + str(i_global), "patch", "none", "false", base)

                append4SymmetricFacestoPatch(patchDict, "fuelTop_" + str(i_global), base1, "top")
        

        ###################################################################################################################
        ####################################  Setting the inner boundary patches ##########################################
        ###################################################################################################################
                
        if rInner>0:
            base=[x + i_vertex for x in [1, 0, shift, shift+1]]
            addToPatchDict(patchDict, "fuelInner", "patch", "none", "false", base)
            addToPatchDict(patchDict, "fuelInner", "patch", "none", "false", [x + 1 for x in base])
            addToPatchDict(patchDict, "fuelInner", "patch", "none", "false", [x + 2 for x in base])
            base=[x + i_vertex for x in [0, 3, shift+3, shift]]
            addToPatchDict(patchDict, "fuelInner", "patch", "none", "false", base)
        

        ###################################################################################################################
        ####################################  Setting the outer boundary patches ##########################################
        ###################################################################################################################

        base=[x + i_vertex for x in [shift-4, shift-3, 2*shift-3, 2*shift-4]]  
        addToPatchDict(patchDict, "fuelOuter", "regionCoupledOFFBEAT", "cladInner", "true", base)
        addToPatchDict(patchDict, "fuelOuter", "regionCoupledOFFBEAT", "cladInner", "true", [x + 1 for x in base])
        addToPatchDict(patchDict, "fuelOuter", "regionCoupledOFFBEAT", "cladInner", "true", [x + 2 for x in base])
        base=[x + i_vertex for x in [shift-1, shift-4, 2*shift-4, 2*shift-1]]
        addToPatchDict(patchDict, "fuelOuter", "regionCoupledOFFBEAT", "cladInner", "true", base)
        



def addCladToPatchDict(patchDict, mergePatchDict, clad_block, merging, N_blocks, i_vertex, i_global, geometry):

    type=clad_block["type"]

    #*********************************************  2D discrete  **********************************************************
    if geometry=="2D-discrete":

        ###################################################################################################################
        ####################################  Setting the bottom boundary patches #########################################
        ###################################################################################################################

        base=[x + i_vertex for x in [0, 1, 3, 2]]
        base1=[x + i_vertex for x in [8, 1, 0, 8]]

        if i_global==1:
            addToPatchDict(patchDict, "cladBottom", "patch", "none", "false", base)
            if(type=="cap"):
                addToPatchDict(patchDict, "cladBottom", "patch", "none", "false", base1)
        else:
            if not(merging):
                addToPatchDict(patchDict, "cladBottom_" + str(i_global), "regionCoupledOFFBEAT", "cladTop_"+ str(i_global-1) , "false", base)
            else:
                addToPatchDict(patchDict, "cladBottom_" + str(i_global), "patch", "none", "false", base)
        
        ###################################################################################################################
        ####################################  Setting the top boundary patches ############################################
        ###################################################################################################################
            
        base=[x + 4 for x in base[::-1]]
        base1=[x + i_vertex for x in [9, 5, 4, 9]]

        if i_global==N_blocks:
            addToPatchDict(patchDict, "cladTop", "patch", "none", "true", base)
            if(type=="cap"):
                addToPatchDict(patchDict, "cladTop", "patch", "none", "true", base1)
        else:
            if not(merging):
                addToPatchDict(patchDict, "cladTop_" + str(i_global), "regionCoupledOFFBEAT", "cladBottom_"+ str(i_global+1) , "false", base)
            else:
                addToPatchDict(patchDict, "cladTop_" + str(i_global), "patch", "none", "true", base)
                mergePatchDict["cladTop_" + str(i_global)] = 'cladBottom_' + str(i_global+1)

        ###################################################################################################################
        ####################################  Setting cap inner patches ###################################################
        ###################################################################################################################

        base1=[x + i_vertex for x in [9, 4, 5, 9]]
        base2=[x + i_vertex for x in [8, 1, 0, 8]]

        if type=="cap":
            if i_global==N_blocks:       
                addToPatchDict(patchDict, "topCapInner", "regionCoupledOFFBEAT", "fuelTop" , "true", base2)
            else:
                addToPatchDict(patchDict, "bottomCapInner", "regionCoupledOFFBEAT", "fuelBottom" , "false", base1)

        ###################################################################################################################
        ####################################  Setting the inner boundary patches ##########################################
        ###################################################################################################################
        if type=="normal":
            base=[x + i_vertex for x in [1, 0, 4, 5]]
            addToPatchDict(patchDict, "cladInner", 'regionCoupledOFFBEAT', 'fuelOuter', "false", base)

        ###################################################################################################################
        ####################################  Setting the outer boundary patches ##########################################
        ###################################################################################################################
        
        base=[x + i_vertex for x in [2, 3, 7, 6]]  
        addToPatchDict(patchDict, "cladOuter", 'patch',  'none', "false", base)

    #**************************************************   3D  **************************************************************       
    else:

        shift=int(clad_block["nVertices"]/2)

        ###################################################################################################################
        ####################################  Setting the bottom boundary patches #########################################
        ###################################################################################################################

        base=[x + i_vertex for x in [0, 3, 2, 1]]
        base1=[x + i_vertex for x in [4, 0, 1, 5]]
        base2=[x + 4 for x in base1]

        if i_global==1:

            addToPatchDict(patchDict, "cladBottom", "patch", "none", "false", [])
            if(type=="cap"):
                addToPatchDict(patchDict, "cladBottom", "patch", "none", "false", base)
                append4SymmetricFacestoPatch(patchDict, "cladBottom", base2, "bottom")

            append4SymmetricFacestoPatch(patchDict, "cladBottom", base1, "bottom")
        else:
            if not(merging):
                addToPatchDict(patchDict, "cladBottom_" + str(i_global), "regionCoupledOFFBEAT", "cladTop_"+ str(i_global-1) , "false", [])
            else:
                addToPatchDict(patchDict, "cladBottom_" + str(i_global), "patch", "none", "false", [])
            
            if(type=="cap"):
                append4SymmetricFacestoPatch(patchDict, "cladBottom_" + str(i_global), base2, "bottom")
            else:
                append4SymmetricFacestoPatch(patchDict, "cladBottom_" + str(i_global), base1, "bottom")

        ###################################################################################################################
        ####################################  Setting the top boundary patches ############################################
        ###################################################################################################################
            
        base=[x + shift for x in base[::-1]]
        base1=[x + shift for x in base1[::-1]]
        base2=[x + shift for x in base2[::-1]]

        if i_global==N_blocks:

            addToPatchDict(patchDict, "cladTop", "patch", "none", "true", [])
            if(type=="cap"):
                addToPatchDict(patchDict, "cladTop", "patch", "none", "true", base)
                append4SymmetricFacestoPatch(patchDict, "cladTop", base2, "top")

            append4SymmetricFacestoPatch(patchDict, "cladTop", base1, "top")
        else:
            if not(merging):
                addToPatchDict(patchDict, "cladTop_" + str(i_global), "regionCoupledOFFBEAT", "cladBottom_"+ str(i_global+1) , "false", [])
            else:
                addToPatchDict(patchDict, "cladTop_" + str(i_global), "patch", "none", "true", [])
                mergePatchDict["cladTop_" + str(i_global)] = 'cladBottom_' + str(i_global+1)

            if type=="cap":
                append4SymmetricFacestoPatch(patchDict, "cladTop_" + str(i_global), base2, "top")
            else:
                append4SymmetricFacestoPatch(patchDict, "cladTop_" + str(i_global), base1, "top")

        ###################################################################################################################
        ####################################  Setting cap inner patches ###################################################
        ###################################################################################################################
        base=[x + i_vertex for x in [0, 3, 2, 1]]
        base1=[x + i_vertex for x in [4, 0, 1, 5]]

        if type=="cap":
            if i_global==N_blocks:       
                    addToPatchDict(patchDict, "topCapInner", "regionCoupledOFFBEAT", "fuelTop" , "true", base)
                    append4SymmetricFacestoPatch(patchDict, "topCapInner", base1, "bottom")
            else:
                    base=[x + shift for x in base[::-1]]
                    base1=[x + shift for x in base1[::-1]]
                    addToPatchDict(patchDict, "bottomCapInner", "regionCoupledOFFBEAT", "fuelBottom" , "false", base)
                    append4SymmetricFacestoPatch(patchDict, "bottomCapInner", base1, "top")

        ###################################################################################################################
        ####################################  Setting the inner boundary patches ##########################################
        ###################################################################################################################
        if type=="normal":
            base=[x + i_vertex for x in [1, 0, shift, shift+1]]
            addToPatchDict(patchDict, "cladInner", 'regionCoupledOFFBEAT', 'fuelOuter', "false", base)
            addToPatchDict(patchDict, "cladInner", 'regionCoupledOFFBEAT', 'fuelOuter', "false", [x + 1 for x in base])
            addToPatchDict(patchDict, "cladInner", 'regionCoupledOFFBEAT', 'fuelOuter', "false", [x + 2 for x in base])
            base=[x + i_vertex for x in [0, 3, shift+3, shift]]
            addToPatchDict(patchDict, "cladInner", 'regionCoupledOFFBEAT', 'fuelOuter', "false", base)
        
        ###################################################################################################################
        ####################################  Setting the outer boundary patches ##########################################
        ###################################################################################################################
        
        base=[x + i_vertex for x in [shift-4, shift-3, 2*shift-3, 2*shift-4]]  
        addToPatchDict(patchDict, "cladOuter", 'patch',  'none', "false", base)
        addToPatchDict(patchDict, "cladOuter", 'patch',  'none', "false", [x + 1 for x in base])
        addToPatchDict(patchDict, "cladOuter", 'patch',  'none', "false", [x + 2 for x in base])
        base=[x + i_vertex for x in [shift-1, shift-4, 2*shift-4, 2*shift-1]]
        addToPatchDict(patchDict, "cladOuter", 'patch',  'none', "false", base)
            

###################################################################################################################################################
#########################----------------------------------- GENERAL WRITING FUNCTIONS -----------------------------------#########################
###################################################################################################################################################
    
def writeHeader(file):
    header = """
/*--------------------------------*- C++ -*----------------------------------*\\
| ========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5.0                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     5.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""
    file.write(header)

def writeGeometry(list_spheres, file):
    if list_spheres:
        file.write("\ngeometry\n{\n")

        for i in range(len(list_spheres)):
            sphere = list_spheres[i]
            file.write(f"\n    sphere_{i}\n")
            file.write("    {\n")
            file.write("        type searchableSphere;\n")
            file.write(f"        centre ({sphere['x']} {sphere['y']} {sphere['z']});\n")
            file.write(f"        radius {sphere['radius']};\n")
            file.write("    }\n")
        file.write("}\n")

def writeVertices(list_vertices, file):
    file.write("\nvertices\n(\n")

    for vertex in list_vertices:
        vertex_str = "    (" + " ".join(map(str, vertex)) + ")\n"
        file.write(vertex_str)

    file.write(");\n")

def writeBlocks(list_blocks, file):
    file.write("\nblocks\n(\n")
    for block in list_blocks:
        vertices_str = " ".join(map(lambda x: str(int(x)), block["vertices"]))
        mesh_str = " ".join(map(lambda x: str(int(x)), block["mesh"]))
        block_str = f"    hex ( {vertices_str} ) {block['name']} ({mesh_str}) simpleGrading (1 1 1)\n"
        file.write(block_str)
    file.write(");\n")

def writeEdges(list_edges, file):
    if list_edges:
        file.write("\nedges\n(\n")
        for edge in list_edges:
            vertices = [int(v) for v in edge['vertices']] # ensuring vertices are treated as integers
            midpoint_str = ' '.join(str(x) for x in edge['midpoint'])
            edge_str = f"    arc {vertices[0]} {vertices[1]} ({midpoint_str})\n"
            file.write(edge_str)
        file.write(");\n")

def writeFaceProjections(list_projection_faces, file):
    if list_projection_faces:
        file.write("\nfaces\n(\n")

        for projection in list_projection_faces:
            face = ' '.join(str(int(x)) for x in projection['face'])
            sphere_name = projection['sphere']
            file.write(f"    project ({face}) {sphere_name}\n")
        file.write(");\n")

def writeBoundaries(patchDict, file):

    file.write("\nboundary\n(\n")
    for patchName, patchInfo in patchDict.items():
        file.write(f"    {patchName}\n")
        file.write("    {\n")
        file.write(f"        type {patchInfo['type']};\n")

        if patchInfo['type'] == "regionCoupledOFFBEAT":
            file.write(f"        neighbourPatch {patchInfo.get('neighbour', '')};\n")
            file.write("        neighbourRegion region0;\n")
            file.write(f"        owner {'true' if patchInfo.get('owner') == 'true' else 'false'};\n")
            # Specific logic for cladInner or fuelOuter
            if patchName == "cladInner" or patchName == "fuelOuter":
                file.write("        updateAMI true;\n")
            else:
                file.write("        updateAMI false;\n")
        
        file.write("        faces\n        (\n")
        for face in patchInfo['faces']:
            faceStr = " ".join(map(str, face))
            file.write(f"            ({faceStr})\n")
        file.write("        );\n")
        file.write("    }\n\n")
    file.write(");\n")

def writeMergedPatches(mergePatchDict, file):
    file.write("\nmergePatchPairs \n(\n")
    for masterPatchName in mergePatchDict:
        slavePatchName = mergePatchDict[masterPatchName]
        file.write("\t(")
        file.write(masterPatchName + " " + slavePatchName)
        file.write(")\n")

    file.write(");\n\n")

'''------------------------------------------------------------
-------------------------- MAIN -------------------------------
------------------------------------------------------------'''
###############################################################
#### Reading the rodDict file and making the dictionary #######
###############################################################
    
# Reading the data from the rodDict file 
with open('rodDict') as f: 
    data = f.read() 

# Reconstructing the data as a dictionary 
rodDict = ast.literal_eval(data)

###############################################################
######### Extracting parameters from the dictionary ###########
###############################################################

# Basic parameters
convertToMeters = rodDict['convertToMeters']
geometry = rodDict['geometryType'] 

nBlocksFuel = rodDict['nBlocksFuel']
blockNameFuel = rodDict['blockNameFuel']
nBlocksClad = rodDict['nBlocksClad']
blockNameClad = rodDict['blockNameClad']

# Geometrical parameters for fuel and cladding
rInnerFuel = rodDict['rInnerFuel']
rOuterFuel = rodDict['rOuterFuel']
rInnerClad = rodDict['rInnerClad']
rOuterClad = rodDict['rOuterClad']
heightFuel = rodDict['heightFuel']
heightClad = rodDict['heightClad']

# Global offsets
offsetFuel = rodDict['offsetFuel']
offsetClad = rodDict['offsetClad']

# Mesh properties for cladding
nCellsZClad = rodDict['nCellsZClad']
nCellsRClad = rodDict['nCellsRClad']

if geometry!='3D':
    # transforming degrees to radians
    wedgeAngle=rodDict['wedgeAngle']*math.pi/180

if geometry=='2D-discrete' or geometry=='3D':
    # Pellet parameters
    nPelletsFuel = rodDict['nPelletsFuel']
    totalPelletNumber=sum(nPelletsFuel)

if geometry=='3D':
    squareFractionTopCap=rodDict.get('squareFractionTopCap', None)
    squareFractionBottomCap=rodDict.get('squareFractionBottomCap', None)

if geometry=='2D-discrete' or geometry=='3D':
    # Global merge patch pairs options
    mergeCladPatchPairs = rodDict['mergeCladPatchPairs']
    mergeFuelPatchPairs = rodDict['mergeFuelPatchPairs']

    # Geometry of the fuel pellet
    rDishFuel = rodDict['rDishFuel']
    rCurvatureDish = rodDict['rCurvatureDish']
    chamferHeight = rodDict['chamferHeight']
    chamferWidth = rodDict['chamferWidth']

if geometry=='3D':
    # Mesh properties for fuel
    squareFraction = rodDict['squareFraction']
    nCellsAzimuthalFuel = rodDict['nCellsAzimuthalFuel']
    nCellsAzimuthalClad = rodDict['nCellsAzimuthalClad']
    eccentricity=rodDict['eccentricity']
    eccentricity_mode=rodDict['eccentricity_mode']

    if eccentricity_mode=='manual':
        eccVector=rodDict['eccentricity_vector']


###################################################################
############## Dealing with the top and bottom caps ###############
####### adding new blocks to the list of cladding blocks ##########
###################################################################
cladType = ['normal' for i in range(nBlocksClad)]

bottomCap=False
topCap=False

if geometry!='1D':
    # Cap parameters
    bottomCapHeight = rodDict['bottomCapHeight']
    topCapHeight = rodDict['topCapHeight']
    nCellsRBottomCap = rodDict.get('nCellsRBottomCap', None)
    nCellsZBottomCap = rodDict.get('nCellsZBottomCap', None)
    nCellsRTopCap = rodDict.get('nCellsRTopCap', None)
    nCellsZTopCap = rodDict.get('nCellsZTopCap', None)
        

    if(bottomCapHeight>0):
        bottomCap=True
        nBlocksClad += 1
        blockNameClad.insert(0, 'cladding')
        cladType.insert(0, 'cap')
        rInnerClad.insert(0, rInnerClad[0])
        rOuterClad.insert(0, rOuterClad[0])
        heightClad.insert(0, bottomCapHeight)
        offsetClad -= float(heightClad[0])
        nCellsRClad.insert(0, nCellsRClad[0])
        nCellsZClad.insert(0, nCellsZBottomCap)
        if geometry=='3D':
            nCellsAzimuthalClad.insert(0, nCellsAzimuthalClad[0])

    if(topCapHeight>0):
        topCap=True
        nBlocksClad += 1
        blockNameClad.append('cladding')
        cladType.append('cap')
        rInnerClad.append(rInnerClad[nBlocksClad-2])
        rOuterClad.append(rOuterClad[nBlocksClad-2])
        heightClad.append(topCapHeight)
        nCellsRClad.append(nCellsRClad[nBlocksClad-2])
        nCellsZClad.append(nCellsZTopCap)
        if geometry== '3D':
            nCellsAzimuthalClad.append(nCellsAzimuthalClad[nBlocksClad-2])

if geometry=='1D' or geometry=='2D-smeared':

    ###########################################################################
    # Initialize lists to hold dictionaries for each fuel and cladding blocks #
    ###########################################################################

    fuel_blocks = []
    cladding_blocks = []

    for i in range(nBlocksFuel):

        block_fuel = {
            
            'name':                      blockNameFuel[i],

            'rInner':                     rInnerFuel[i],
            'rOuter':                    rOuterFuel[i],
            'height':                    heightFuel[i],

            # mesh properties:
            "nR":                        rodDict['nCellFuelR'][i],
            "nZ":                        rodDict['nCellFuelZ'][i],

            "nVertices":                 8

        }

        fuel_blocks.append(block_fuel)


    for i in range(nBlocksClad):
        block_clad = {

            'name':                       blockNameClad[i],
            'rInner':                     rInnerClad[i],
            'rOuter':                     rOuterClad[i],
            'height':                     heightClad[i],

            'type':                       cladType[i],

            "nR":                         nCellsRClad[i],
            "nZ":                         nCellsZClad[i],

            "nVertices":                 8
        }    

        if cladType[i]=='cap':
            block_clad["nVertices"]=10
            if i==0:
                block_clad["nRInner"]=nCellsRBottomCap
            else:
                block_clad["nRInner"]=nCellsRTopCap

        cladding_blocks.append(block_clad)
    

if geometry=='2D-discrete' or geometry=='3D':

    nCellsRPellet=rodDict.get('nCellsRPellet',None)
    nCellsRDish = rodDict.get('nCellsRDish', None)
    nCellsRChamfer = rodDict.get('nCellsRChamfer',None)
    nCellsZPellet = rodDict['nCellsZPellet']

    # Creating the list of pellet types
    pelletType = []
    nVerticesFuel = [0.0 for i in range(nBlocksFuel)]
    rLandFuel     = [0.0 for i in range(nBlocksFuel)] # end of land / start of chamfer

    for i in range(nBlocksFuel):
        if rDishFuel[i] > 0.0 and chamferWidth[i] > 0.0:
            pelletType.append('dishedChamfered')
            rLandFuel[i]     = rOuterFuel[i] - chamferWidth[i]
            if geometry=='3D':
                nVerticesFuel[i] = 32
            else:
                nVerticesFuel[i] = 16
        elif rDishFuel[i] > 0.0:
            pelletType.append('dished')
            rLandFuel[i]     = rOuterFuel[i]
            if geometry=='3D':
                nVerticesFuel[i] = 24
            else:
                nVerticesFuel[i] = 12
        elif chamferWidth[i] > 0.0:
            pelletType.append('chamfered')
            rLandFuel[i]     = rOuterFuel[i] - chamferWidth[i]
            if geometry=='3D':
                nVerticesFuel[i] = 24
            else:
                nVerticesFuel[i] = 12
        else:
            pelletType.append('flat')
            rLandFuel[i]     = rOuterFuel[i]
            if geometry=='3D':
                nVerticesFuel[i] = 16
            else:
                nVerticesFuel[i]= 8


    '''------------------------------------------------------------------------
    ------------------------------------------------------------------------'''
    ###########################################################################
    # Initialize lists to hold dictionaries for each fuel and cladding blocks #
    ###########################################################################
    
    fuel_blocks = []
    cladding_blocks = []

    for i in range(nBlocksFuel):
        fuel_block = {

        "blockName": blockNameFuel[i],

        "rInner": rInnerFuel[i],
        "rOuter": rOuterFuel[i],
        "height": heightFuel[i]/nPelletsFuel[i], # height of each pellet

        "type": pelletType[i],

        "rLand": rLandFuel[i],
        "rDish": rDishFuel[i],
        "rCurvatureDish" : rCurvatureDish[i],
        "chamferWidth": chamferWidth[i],
        "chamferHeight": chamferHeight[i],

        "nCellsRPellet": nCellsRPellet[i],
        "nCellsRDish": nCellsRDish[i],
        "nCellsRChamfer": nCellsRChamfer[i],    
        "nCellsZPellet": nCellsZPellet[i],
        "nVertices": nVerticesFuel[i],
        }

        if geometry=='3D':
            fuel_block['squareFraction'] = squareFraction[i]
            fuel_block["nCellsAzimuthal"]= nCellsAzimuthalFuel[i]
        else:
            fuel_block['wedgeAngle'] = wedgeAngle
      
        fuel_blocks.append(fuel_block)

    for i in range(nBlocksClad):
        clad_block = {

            "blockName": blockNameClad[i],

            "type": cladType[i],
            "rInner": rInnerClad[i],
            "rOuter": rOuterClad[i],
            "height": heightClad[i],
            "nVertices": 8,
            "nCellsR": nCellsRClad[i],
            "nCellsZ": nCellsZClad[i]
        }

        if geometry=='3D':
            clad_block["nCellsAzimuthal"]= nCellsAzimuthalClad[i]
            clad_block["nVertices"]= 16
        else:
            clad_block['wedgeAngle'] = wedgeAngle

        if cladType[i]=="cap":
            if i==0:
                clad_block["nCellsRInner"]=nCellsRBottomCap
                if geometry=='3D':
                    clad_block["squareFraction"]=squareFractionBottomCap
            else:
                clad_block["nCellsRInner"]=nCellsRTopCap
                if geometry=='3D':
                    clad_block["squareFraction"]=squareFractionTopCap
            if geometry=='3D':
                clad_block["nVertices"]=24
            else:
                clad_block["nVertices"]=10

        cladding_blocks.append(clad_block)


global_clad_offset=offsetClad
global_fuel_offset=offsetFuel

list_spheres=[]
list_vertices=[]
list_blocks=[]
list_projection_faces=[]
list_edges=[]

patchDict = {}
mergePatchDict = defaultdict(list)

i_global=1
i_vertex=0
i_sphere=2


if geometry=="1D" or geometry=="2D-smeared":
        
        for i in range(nBlocksFuel):
            addWedgeVertices(list_vertices, fuel_blocks[i], wedgeAngle, global_fuel_offset, 0)
            addWedgeBlocks(list_blocks, fuel_blocks[i], i_vertex, 0)
            addWedgePatches(patchDict, mergePatchDict, fuel_blocks[i], nBlocksFuel, bottomCap, topCap, i_vertex, i_global, geometry, 0)  
            global_fuel_offset+=fuel_blocks[i]['height']
            i_vertex+=fuel_blocks[i]['nVertices']
            i_global+=1

        i_global=1

        for i in range(nBlocksClad):
            addWedgeVertices(list_vertices, cladding_blocks[i], wedgeAngle, global_clad_offset, 1)
            addWedgeBlocks(list_blocks, cladding_blocks[i], i_vertex, 1)
            addWedgePatches(patchDict, mergePatchDict, cladding_blocks[i], nBlocksClad, bottomCap, topCap, i_vertex, i_global, geometry, 1)
            global_clad_offset+=cladding_blocks[i]['height']
            i_vertex+=cladding_blocks[i]['nVertices']
            i_global+=1
            



# finding the minimum gap width along the whole rod, used for the case when
# simulating randoom pellet eccentricity in 3D case by 'default model'
if geometry=='3D' and eccentricity and eccentricity_mode=='default':
    max_rOuterFuel=0
    for i in range(nBlocksFuel):
        rOuter=fuel_blocks[i]["rOuter"]
        if rOuter>max_rOuterFuel:
            max_rOuterFuel=rOuter

    min_RInnerClad=100000
    for i in range(nBlocksClad):
        rInner=cladding_blocks[i]["rInner"]
        if rInner<min_RInnerClad:
            min_RInnerClad=rInner

    minGap=min_RInnerClad-max_rOuterFuel

# Call the function with your desired file name
if geometry=="3D" or geometry=="2D-discrete":

    shiftX=0
    shiftY=0

    for i in range(nBlocksFuel):
        for j in range(nPelletsFuel[i]):
           
            if geometry=='3D': 
                if eccentricity:
                    if eccentricity_mode=='default':
                        shiftX=random.uniform(-minGap, minGap)
                        shiftY=random.uniform(-minGap, minGap)
                    else:
                        shiftX=eccVector[i_global-1][0]
                        shiftY=eccVector[i_global-1][1]             

            addSpheres(list_spheres,fuel_blocks[i], global_fuel_offset, geometry, shiftX, shiftY)
            addPelletVertices(list_vertices, fuel_blocks[i], global_fuel_offset, geometry, shiftX, shiftY)
            addFuelBlocks(list_blocks, fuel_blocks[i], i_vertex, geometry)
            addFuelEdges(list_edges, fuel_blocks[i], i_vertex, global_fuel_offset, geometry, shiftX, shiftY)
            i_sphere = addFaceProjections(list_projection_faces, fuel_blocks[i], i_vertex, i_sphere, geometry)
            addFuelToPatchDict(patchDict, mergePatchDict, fuel_blocks[i], mergeFuelPatchPairs, totalPelletNumber, bottomCap, topCap, i_vertex, i_global, geometry)
            i_vertex+=fuel_blocks[i]['nVertices']
            global_fuel_offset+=fuel_blocks[i]['height']
            i_global+=1


    i_global=1
    for i in range(nBlocksClad):
            addCladVertices(list_vertices, cladding_blocks[i], global_clad_offset, geometry)
            addCladBlocks(list_blocks, cladding_blocks[i], i_vertex, geometry)
            addCladEdges(list_edges, cladding_blocks[i], i_vertex, global_clad_offset, geometry)
            addCladToPatchDict(patchDict, mergePatchDict, cladding_blocks[i], mergeCladPatchPairs, nBlocksClad, i_vertex, i_global, geometry)
            i_vertex+=cladding_blocks[i]['nVertices']
            global_clad_offset+=cladding_blocks[i]['height']
            i_global+=1




##################################################################
##### Now, all od the parameters are set up...####################
## The only thing left is to write them in the blockMeshDict :) ##
##################################################################

file = open("blockMeshDict", "w+")
writeHeader(file)
file.write("\nconvertToMeters " + str(convertToMeters) + "; \n\n")
writeGeometry(list_spheres,file)
writeVertices(list_vertices,file)
writeBlocks(list_blocks, file)
writeEdges(list_edges, file)
writeFaceProjections(list_projection_faces, file)
writeBoundaries(patchDict, file)
writeMergedPatches(mergePatchDict, file)

'''
# Function to print variables to a file
def print_variables_to_file(filename):
    with open(filename, 'w') as file:
        for name, value in globals().items():
            # Optionally, filter out built-in variables and modules
            if not name.startswith('__') and not hasattr(value, '__call__'):
                file.write(f"{name}: {value}\n")

# Call the function with your desired file name
print_variables_to_file('variables.txt')
'''