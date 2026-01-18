# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Q-SLBL
                                 A QGIS plugin
 Computes the estimated landslide volume using SLBL
                             -------------------
        begin                : 2026-01-12
        copyright            : (C) 2026 by Richmond Camingawan
        email                : rdtcamingawan@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""

def classFactory(iface):
    # load StationLines class from file StationLines
    from .stationlines import StationLinesPlugin
    return QSLBL(iface)
