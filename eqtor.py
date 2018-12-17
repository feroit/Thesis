from math import *
import numpy as np
import bpy
import time
import ctypes
import _ctypes
class param:
    def __init__(self, sizet, sizex, sizey, R, r, vmax, vmin, dv, time, a, name, u0, f):
        self.sizet = sizet
        self.sizex = sizex
        self.sizey = sizey
        self.R = R
        self.r = r
        self.vmax = vmax
        self.vmin = vmin
        self.dv = dv
        self.time = time
        self.a = a
        self.name = name
        self.data = np.zeros((sizet, sizex*sizey), dtype=np.double)
        self.u0 = u0
        self.f = f
    def __del__(self):
        self.data = None
    def KtoRGB(self, v):#Функция перевода температуры в цвет
        vmin = self.vmin
        vmax = self.vmax
        dv = self.dv
        r, g, b = 1., 1., 1.
        if v < vmin:
            v = vmin
        if v > vmax:
            v = vmax
        if v < (vmin + 0.25 * dv):
            r = 0.
            g = 4. * (v - vmin) / dv
        elif v < (vmin + 0.5 * dv):
            r = 0.
            b = 1. + 4. * (vmin + 0.25 * dv - v) / dv
        elif v < (vmin + 0.75 * dv):
            r = 4. * (v - vmin - 0.5 * dv) / dv
            b = 0.
        else:
            g = 1. + 4. * (vmin + 0.75 * dv - v) / dv
            b = 0.
        return [r, g, b]
    def KtoRadius(self,v):#Функция перевода температуры в радиус
        if v < self.vmin:
            v = self.vmin
        if v > self.vmax:
            v = self.vmax
        v = (v - self.vmin)/self.dv
        v = 0.5 + v*(self.R-1.)
        return v
    def changecolor(self, vertex_colors, u):#Изменение цвета
        color_map_collection = vertex_colors
        if len(color_map_collection)==0:
            color_map_collection.new()
        color_map = color_map_collection['Col']
        sizex = self.sizex
        sizey = self.sizey
        k = 0
        h_range = range(1,sizex)
        w_range = range(1,sizey)
        for i in h_range:
            for j in w_range:
                one = j-1+(i-1)*sizey
                two = j+(i-1)*sizey
                thr = j+i*sizey
                fou = j-1+i*sizey
                color_map.data[k].color = self.KtoRGB(float(u[one]))
                color_map.data[k+1].color = self.KtoRGB(float(u[two]))
                color_map.data[k+2].color = self.KtoRGB(float(u[thr]))
                color_map.data[k+3].color = self.KtoRGB(float(u[fou]))
                k +=4
    def changecoo(self,vertices, u):#Изменение координат
        sizex = self.sizex
        sizey = self.sizey
        R = self.R
        h_range = range(0,sizex)
        w_range = range(0,sizey)
        dn = (2.0*pi)/(sizex-1.0)
        dm = (2.0*pi)/(sizey-1.0)
        for i in h_range:
            for j in w_range:
                _r = self.KtoRadius(u[j+i*sizey])
                vertices[j+i*sizey].co = (R+_r*cos(i*dn))*cos(j*dm), (R+_r*cos(i*dn))*sin(j*dm), _r*sin(i*dn)
    class Profiler(object):#Класс для замеров времени работы частей кода
        def __enter__(self):
            self._startTime = time.time()
         
        def __exit__(self, type, value, traceback):
            print("Elapsed time: {:.3f} sec".format(time.time() - self._startTime))
    def Recalc(self, value, coobool, colorbool):# Функция изменения параметров меша в зависимости от кадра
        u = self.data[value]
        sizex = self.sizex
        sizey = self.sizey
        mesh = bpy.data.objects['Poly'].data
        if coobool:
            self.changecoo(mesh.vertices, u)
        if colorbool:
            self.changecolor(mesh.vertex_colors, u)
    def Solution(self):# Функция взаимодействия с динамической библиотекой для нахождения решения
        data = self.data
        sizet = self.sizet
        sizex = self.sizex
        sizey = self.sizey
        u0 = self.u0
        f = self.f
        arr2 = [sizet, sizex - 1, sizey - 1]
        arr =[self.R, self.r, self.time, self.a]
        ctypes_arrays = [np.ctypeslib.as_ctypes(array) for array in data]
        point = (ctypes.c_double * 4)(*arr)
        point2 = (ctypes.c_int * 3)(*arr2)
        pointer_ar = (ctypes.POINTER(ctypes.c_double) * sizet)(*ctypes_arrays)
        callback_type = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double)
        callback_type2 = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
        mydll = ctypes.CDLL("./DipDll.dll")
        mydll.Solution(pointer_ar, callback_type(u0), callback_type2(f), point, point2)
        _ctypes.FreeLibrary(mydll._handle)
        del mydll
    def getvert(self, R, r, sizex, sizey):#Создание массива точек тора
        dphi = 2.0*pi/(sizex-1)
        dpsi = 2.0*pi/(sizey-1)
        v = []
        for i in range(0, sizex):
            for j in range(0, sizey):
                v.append([(R+r*cos(dphi*i))*cos(dpsi*j), (R+r*cos(dphi*i))*sin(dpsi*j), r*sin(dphi*i)])
        return v
    def getfac(self, sizex,sizey):#Функция создание связей сетки меша
        x = range(1, sizex)
        y = range(1, sizey)
        fac = []
        for xi in x:
            for yi in y:
                one = yi-1+(xi-1)*sizey
                two = yi+(xi-1)*sizey
                thr = yi+xi*sizey
                fou = yi-1+xi*sizey
                fac.append((one, two, thr, fou))
        return fac
    def Create(self):#Создание тора и его материала
        sizex = self.sizex
        sizey = self.sizey
        R = self.R
        r = self.r
        name = self.name
        fac = self.getfac(sizex, sizey)
        vert = self.getvert(R, r, sizex, sizey)
        f_change = bpy.app.handlers.frame_change_pre
        del f_change[0:len(f_change)]
        ob1 = bpy.data.objects
        if name in ob1:
            ob1.remove(ob1[name], True)
        mesh = bpy.data.meshes.new('Mesh')
        ob = bpy.data.objects.new(name, mesh)
        ob.location = (0,0,0)
        ob.show_name = True
        bpy.context.scene.objects.link(ob)
        mesh.from_pydata(vert, [], fac)
        mesh.update(calc_edges=True)
        mat = bpy.data.materials.new('vertex_material')
        mat.use_vertex_color_paint = True
        mat.use_vertex_color_light = True  # material affected by lights
        ob.data.materials.append(mat)