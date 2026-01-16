"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

from typing import Any, Optional

from qgis.core import (
    QgsFeatureSink,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingContext,
    QgsProcessingException,
    QgsProcessingFeedback,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterEnum,
    QgsProcessingParameterNumber,
)
from qgis import processing

# Processing imports
import numpy as np
import deepcopy


class ExampleProcessingAlgorithm(QgsProcessingAlgorithm):
    """
    This is an example algorithm that takes a vector layer and
    creates a new identical one.

    It is meant to be used as an example of how to create your own
    algorithms and explain methods and variables used to do it. An
    algorithm like this will be available in all elements, and there
    is not need for additional work.

    All Processing algorithms should extend the QgsProcessingAlgorithm
    class.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    def SLBL(self,grid_dem,grid_mask,tol,maxt,maxv,z_min,planes=None):
	# initilizes thickness and difference grids
	s=grid_dem.shape
	grid_thickn = np.zeros(s) # This just creates an array of zer0s with the shape equal to the grid_dem
	grid_diff = np.ones(s) # This just creates an array of ones (1) with the shape equal to the grid_dem

	# Creates a matrice to store the values of the neighbouring cells in the previous iteration
	mat_neigh = np.zeros(s)
	mat_neigh = np.expand_dims(mat_neigh,axis=2) # just changes the dimensionality of the array. no changes in the data, just the shape of the array. 

	# NOTE: nb_neigh is initialized at the main script
	# nb_neigh is 4, if the method selected uses "4 neighbors, average" or "4 neighbors, min/max" 
	if nb_neigh ==4:
		mat_neigh = np.tile(mat_neigh,(1,1,4)) # np.tile repeats plane x1, rows x1, columns x4. This creates a 4 copies of the 
	else: 
		mat_neigh = np.tile(mat_neigh,(1,1,8))

	# Creates a matrice where the proposed value and previous value are stored for comparison
	mat_comp = np.zeros(s)
	mat_comp = np.expand_dims(mat_comp,axis=2)
	mat_comp = np.tile(mat_comp,(1,1,2))

	# Initiate the slbl grid (altitude)
	grid_slbl = deepcopy(grid_dem) # RDTC: creates a copy and make sure nothing happens to the original

	nb_iter = 0
	volume = 0.

	## RDTC:
	# Tests if maxt (max thickness) is a finite number.
	# If True, then the original DEM is substracted by the MAXT
	if np.isfinite(maxt):
		grid_maxt = grid_dem - maxt

	# The SLBL strarts here
	while np.amax(grid_diff)>stop and volume<maxv:
		nb_iter=nb_iter+1
		grid_thickn_prev = deepcopy(grid_thickn) # rdtc: just an array of 0's at the start
		grid_slbl_prev = deepcopy(grid_slbl) # rdtc: the original DEM at the start
		
		# writes the values of the neighbourings cells in the 3rd dimension of the matrix
		mat_neigh[:-1,:,0]=grid_slbl_prev[1:,:] # South Neighbor
		mat_neigh[1:,:,1]=grid_slbl_prev[:-1,:] # North Neighbor
		mat_neigh[:,:-1,2]=grid_slbl_prev[:,1:] # East Neighbor
		mat_neigh[:,1:,3]=grid_slbl_prev[:,:-1] # West Neighbor
		
		# diagonals
		# rdtc: only acessed when method is 8-neighbors
		if nb_neigh ==8:
			mat_neigh[:-1,:-1,4]=grid_slbl_prev[1:,1:] # South East Neighbor
			mat_neigh[:-1,1:,5]=grid_slbl_prev[1:,:-1] # South West Neighbor
			mat_neigh[1:,1:,6]=grid_slbl_prev[:-1,:-1] # North West Neighbor
			mat_neigh[1:,:-1,7]=grid_slbl_prev[:-1,1:] # North East Neigbor
		
		# rdtc: default criteria is minmax, when method is 4-neighbors, minmax
		# rdtc: 4N can also use average, as well as, 8N
		if criteria == 'minmax':
			mat_max=np.amax(mat_neigh,axis=2)
			mat_min=np.amin(mat_neigh,axis=2)

			# rdtc: mat_mean is the z_temp in Jaboyedoff's
			mat_mean=(mat_max+mat_min)/2

		elif criteria == 'average':
			mat_mean=np.mean(mat_neigh,axis=2)

		# rdtc: add tolerance for a more parabolic slip surface
		mat_mean=mat_mean+tol

		# limit to the lower altitude around the polygon
		# rdtc: checks whether current mean is lower than the defined Z_min
		# if yes, then Z_min would be used
		if np.isfinite(z_min):
			mat_mean=np.maximum(mat_mean,z_min)

		# limit to the maximum thickness
		# same logic as Z_min
		if np.isfinite(maxt):
			mat_mean=np.maximum(mat_mean,grid_maxt)
		
		# rdtc: only accessed when z_min and maxt are infiinte numbers
		else not planes is None:
			mat_mean=np.maximum(mat_mean,planes)
		
		# rdtc:Updates the values of mat_comp matrix
		# plane 0: values of mat mean
		# plane 1: values of grid_slbl
		mat_comp[:,:,0]=mat_mean
		mat_comp[:,:,1]=grid_slbl
		
		# Check if the new value should be kept
		# rdtc: what is this inverse? 
		if inverse == 'true':
			grid_slbl=np.amax(mat_comp,axis=2)
		else:
			grid_slbl=np.amin(mat_comp,axis=2)
		
		# Replaces the values of the SLBL by the original values outside the masked area
		grid_slbl[~grid_mask]=grid_dem[~grid_mask]
		
		# rdtc: This is the difference between the original dem - updated SLBL surface
		grid_thickn = np.absolute(grid_dem - grid_slbl)

		# rdtc: This is the difference between the previous and current iteration
		# of the grid
		grid_diff = np.absolute(grid_thickn - grid_thickn_prev)
		
		# rdtc: computation of volume
		volume = (np.sum(grid_thickn)*cellSize*cellSize)
		
		if nb_iter%100==0:
			str_message = '{} iterations. Max diff is {:.3e} m, max thickness is {:.3f} m and volume is {:.6f} million m3'.format(nb_iter,np.amax(grid_diff),np.amax(grid_thickn),volume/1000000)
			QgsProcessingFeedback.pushInfo(str_message)
	# The SLBL is finished
	
	return grid_slbl, grid_thickn, nb_iter

    DEM = "INPUT"
    MASKING = 'MASKING'
    TOL_MODE = 'TOLERANCE_MODE'
    MAXT = 'MAXT'
    MAXV = 'MAXV'
    METHOD = 'METHOD'
    OUTPUT = "OUTPUT"

    def initAlgorithm(self, config: Optional[dict[str, Any]] = None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        # Input required layers
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.DEM,
                "Input DEM Layer",
                [QgsProcessing.SourceType.TypeRaster],
            )
        )
        
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.MASKING,
                "Input Masking Layer",
                [QgsProcessing.SourceType.TypeVectorAnyGeometry],
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.TOL_MODE,
                "Input Tolerance Mode",
                options = [
                    'Single Value',
                    'Auto',
                    'Auto min / inter / max'
                ],
                defaultValue = 0,
                optional=False            
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAXT,
                'Input Maximum Thickness / Depth',
                type = QgsProcessingParameterNumber.Double,
                minValue=0,
                optional = False
        ))
        
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAXV,
                'Input Maximum Volume',
                type = QgsProcessingParameterNumber.Double,
                minValue=0,
                optional=False
            ))
        
        self.addParameter(
            QgsProcessingParameterEnum(
                self.METHOD,
                "Input Method Type",
                options = [
                    '4 neighbours, min/max',
                    '8 neighbours, min/max',
                    '4 neighbours, average',
                    '8 neighbours, average',
                    ],
                defaultValue = '4 neighbours, min/max',
            )
        )
        # We add a feature sink in which to store our processed features (this
        # usually takes the form of a newly created vector layer when the
        # algorithm is run in QGIS).
        self.addParameter(
            QgsProcessingParameterFeatureSink(self.OUTPUT, "Output layer")
        )

    def processAlgorithm(
        self,
        parameters: dict[str, Any],
        context: QgsProcessingContext,
        feedback: QgsProcessingFeedback,
    ) -> dict[str, Any]:
        """
        Here is where the processing itself takes place.
        """

        # Retrieve the feature source and sink. The 'dest_id' variable is used
        # to uniquely identify the feature sink, and must be included in the
        # dictionary returned by the processAlgorithm function.
        source = self.parameterAsSource(parameters, self.INPUT, context)

        # If source was not found, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSourceError method to return a standard
        # helper text for when a source cannot be evaluated
        if source is None:
            raise QgsProcessingException(
                self.invalidSourceError(parameters, self.INPUT)
            )

        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            source.fields(),
            source.wkbType(),
            source.sourceCrs(),
        )

        # Send some information to the user
        feedback.pushInfo(f"CRS is {source.sourceCrs().authid()}")

        # If sink was not created, throw an exception to indicate that the algorithm
        # encountered a fatal error. The exception text can be any string, but in this
        # case we use the pre-built invalidSinkError method to return a standard
        # helper text for when a sink cannot be evaluated
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))

        # Compute the number of steps to display within the progress bar and
        # get features from source
        total = 100.0 / source.featureCount() if source.featureCount() else 0
        features = source.getFeatures()

        for current, feature in enumerate(features):
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break

            # Add a feature in the sink
            sink.addFeature(feature, QgsFeatureSink.Flag.FastInsert)

            # Update the progress bar
            feedback.setProgress(int(current * total))

        # To run another Processing algorithm as part of this algorithm, you can use
        # processing.run(...). Make sure you pass the current context and feedback
        # to processing.run to ensure that all temporary layer outputs are available
        # to the executed algorithm, and that the executed algorithm can send feedback
        # reports to the user (and correctly handle cancellation and progress reports!)
        if False:
            buffered_layer = processing.run(
                "native:buffer",
                {
                    "INPUT": dest_id,
                    "DISTANCE": 1.5,
                    "SEGMENTS": 5,
                    "END_CAP_STYLE": 0,
                    "JOIN_STYLE": 0,
                    "MITER_LIMIT": 2,
                    "DISSOLVE": False,
                    "OUTPUT": "memory:",
                },
                context=context,
                feedback=feedback,
            )["OUTPUT"]

        # Return the results of the algorithm. In this case our only result is
        # the feature sink which contains the processed features, but some
        # algorithms may return multiple feature sinks, calculated numeric
        # statistics, etc. These should all be included in the returned
        # dictionary, with keys matching the feature corresponding parameter
        # or output names.
        return {self.OUTPUT: dest_id}
	
    def name(self) -> str:
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return "myscript"

    def displayName(self) -> str:
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return "My Script"

    def group(self) -> str:
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return "Example scripts"

    def groupId(self) -> str:
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return "examplescripts"

    def shortHelpString(self) -> str:
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it.
        """
        return "Example algorithm short description"

    def createInstance(self):
        return self.__class__()
