from paraview.simple import *

reader = LegacyVTKReader(FileNames=["K_full.vtk"])
Show(reader)
rv = GetActiveViewOrCreate("RenderView")
ResetCamera()

# -------------------------------------------------------------
# K_0 : slice + color map + opacité
# -------------------------------------------------------------
calcK0 = Calculator(Input=reader)
calcK0.ResultArrayName = "K_0_diagonal"
calcK0.Function = "K_0"

sliceK0 = Slice(Input=calcK0)
sliceK0.SliceType = "Plane"
sliceK0.SliceType.Origin = [0, 0, 0]
sliceK0.SliceType.Normal = [0, 1, 0]
Show(sliceK0, rv)
sliceK0Rep = GetDisplayProperties(sliceK0, view=rv)
sliceK0Rep.Representation = "Surface"
ColorBy(sliceK0Rep, ('POINTS', 'K_0_diagonal'))
sliceK0Rep.SetScalarBarVisibility(rv, True)
sliceK0Rep.Opacity = 0.3

k0_ctf = GetColorTransferFunction('K_0_diagonal')
k0_ctf.ApplyPreset('Cool to Warm (Extended)', True)
k0_ctf.RescaleTransferFunction(-0.2, 0.2)

# -------------------------------------------------------------
# K_4 : volume + opacité "peak" + autre palette
# -------------------------------------------------------------
calcK4 = Calculator(Input=reader)
calcK4.ResultArrayName = "K_4_diagonal"
calcK4.Function = "K_4"
Show(calcK4, rv)
calcK4Rep = GetDisplayProperties(calcK4, view=rv)
calcK4Rep.Representation = "Volume"
ColorBy(calcK4Rep, ('POINTS', 'K_4_diagonal'))
calcK4Rep.SetScalarBarVisibility(rv, True)

k4_ctf = GetColorTransferFunction('K_4_diagonal')
k4_ctf.ApplyPreset('RdBu', True)
k4_ctf.RescaleTransferFunction(-0.2, 0.2)

k4_otf = GetOpacityTransferFunction('K_4_diagonal')
k4_otf.Points = [
    -0.2, 0.0, 0.5, 0.0,
    -0.05, 0.7, 0.5, 0.0,
    0.0, 1.0, 0.5, 0.0,
    0.05, 0.7, 0.5, 0.0,
    0.2, 0.0, 0.5, 0.0
]

# -------------------------------------------------------------
# K_8 : isosurfaces + color map
# -------------------------------------------------------------
calcK8 = Calculator(Input=reader)
calcK8.ResultArrayName = "K_8_diagonal"
calcK8.Function = "K_8"
contK8 = Contour(Input=calcK8)
contK8.ContourBy = ['POINTS', 'K_8_diagonal']
contK8.Isosurfaces = [-0.1, 0.0, 0.1, 0.2]
Show(contK8, rv)
contK8Rep = GetDisplayProperties(contK8, view=rv)
ColorBy(contK8Rep, ('POINTS', 'K_8_diagonal'))
contK8Rep.SetScalarBarVisibility(rv, True)
contK8Rep.Opacity = 0.3

k8_ctf = GetColorTransferFunction('K_8_diagonal')
k8_ctf.ApplyPreset('Rainbow Desaturated', True)
k8_ctf.RescaleTransferFunction(-0.2, 0.2)

# -------------------------------------------------------------
# Trace : K_0 + K_4 + K_8 en volume
# -------------------------------------------------------------
calcTrace = Calculator(Input=reader)
calcTrace.ResultArrayName = "K_trace"
calcTrace.Function = "K_0 + K_4 + K_8"
Show(calcTrace, rv)
traceRep = GetDisplayProperties(calcTrace, view=rv)
traceRep.Representation = "Volume"
ColorBy(traceRep, ('POINTS', 'K_trace'))
traceRep.SetScalarBarVisibility(rv, True)
traceRep.Opacity = 0.3

trace_ctf = GetColorTransferFunction('K_trace')
trace_ctf.ApplyPreset('Viridis (matplotlib)', True)
trace_ctf.RescaleTransferFunction(-0.3, 0.3)

trace_otf = GetOpacityTransferFunction('K_trace')
trace_otf.Points = [
    -0.3, 0.0, 0.5, 0.0,
    -0.1, 0.6, 0.5, 0.0,
    0.0, 1.0, 0.5, 0.0,
    0.1, 0.6, 0.5, 0.0,
    0.3, 0.0, 0.5, 0.0
]

ResetCamera()
Render()
Interact()
