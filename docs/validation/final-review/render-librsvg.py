import ctypes as C
from pathlib import Path
import sys
fontconfig=C.CDLL('libfontconfig.so.1')
fontconfig.FcConfigGetCurrent.restype=C.c_void_p
fontconfig.FcConfigAppFontAddFile.argtypes=[C.c_void_p,C.c_char_p]
fontconfig.FcConfigAppFontAddFile.restype=C.c_int
fontconfig.FcConfigBuildFonts.argtypes=[C.c_void_p]; fontconfig.FcConfigBuildFonts.restype=C.c_int
cfg=fontconfig.FcConfigGetCurrent()
assert fontconfig.FcConfigAppFontAddFile(cfg,str(Path(sys.argv[3]).resolve()).encode())
assert fontconfig.FcConfigBuildFonts(cfg)
cairo=C.CDLL('libcairo.so.2'); svg=C.CDLL('librsvg-2.so.2'); gobj=C.CDLL('libgobject-2.0.so.0')
cairo.cairo_image_surface_create.argtypes=[C.c_int,C.c_int,C.c_int]; cairo.cairo_image_surface_create.restype=C.c_void_p
cairo.cairo_create.argtypes=[C.c_void_p]; cairo.cairo_create.restype=C.c_void_p
svg.rsvg_handle_new_from_file.argtypes=[C.c_char_p,C.POINTER(C.c_void_p)]; svg.rsvg_handle_new_from_file.restype=C.c_void_p
svg.rsvg_handle_render_cairo.argtypes=[C.c_void_p,C.c_void_p]; svg.rsvg_handle_render_cairo.restype=C.c_int
cairo.cairo_surface_write_to_png.argtypes=[C.c_void_p,C.c_char_p]; cairo.cairo_surface_write_to_png.restype=C.c_int
cairo.cairo_destroy.argtypes=[C.c_void_p]; cairo.cairo_surface_destroy.argtypes=[C.c_void_p]; gobj.g_object_unref.argtypes=[C.c_void_p]
err=C.c_void_p(); handle=svg.rsvg_handle_new_from_file(str(Path(sys.argv[1]).resolve()).encode(),C.byref(err)); assert handle
surface=cairo.cairo_image_surface_create(0,3000,3000); cr=cairo.cairo_create(surface)
assert svg.rsvg_handle_render_cairo(handle,cr)
assert cairo.cairo_surface_write_to_png(surface,str(Path(sys.argv[2]).resolve()).encode())==0
cairo.cairo_destroy(cr); cairo.cairo_surface_destroy(surface); gobj.g_object_unref(handle)
print('librsvg API rendered with the candidate font registered in this process')
