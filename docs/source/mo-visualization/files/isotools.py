from typing import Optional

def append_isosurface(
        pipeline: Pipeline, 
        isolevel: Optional[float] = 0.1,
        color: Optional[tuple[float,float,float]] = None, 
        alpha: Optional[float] = 0.0) -> Pipeline:
    """Append a modifier to create an isosurface from volumetric data
    loaded from a FHI-aims cube file.

    Args:
        pipeline: an OVITO pipeline containing the volumetric data.
        isolevel: numerical value at which to draw the isosruface. Defaults to 0.1.
        color: isosruface color specified as an RGB tuple. Defaults to red.
        alpha: transparency of the isosurface. Defaults to 0.0.

    Returns:
        Pipeline: original pipeline with the modifier appended to it.
    """
    from ovito.vis import SurfaceMeshVis
    from ovito.modifiers import CreateIsosurfaceModifier
    # Create an object to control the visualization of the surface
    smv = SurfaceMeshVis()
    # Set the surface to be uniformly coloured
    smv.color_mapping_mode = SurfaceMeshVis.ColorMappingMode.Uniform
    if color is None:
        color = (1.0, 0.0, 0.0) # red
    # colour value as RGB tuple
    smv.surface_color = color
    # transparency level
    smv.surface_transparency = alpha
    if isolevel < 0.0:
        # invert direction of surface for negative values 
        # otherwise everything from the surface *outwards* gets shaded
        smv.reverse_orientation = True
    pipeline.modifiers.append(CreateIsosurfaceModifier(
        operate_on = "voxels:imported",
        isolevel = isolevel,
        property = 'Property',
        vis = smv))
    return pipeline
