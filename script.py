import sys
sys.path.append("./")
from eqtor import *
_frameold = 0
T, N, M = 720, 100, 200 #Сетка по времени, малый и большой круг
R, r = 8., 4. #Малый и большой радиус
time, a = 10., 1.0 #размер временного интервала, коэффициент температуропроводности
tmax, tmin, dt = 9., -9., 18. #Параметры отрисовки
p = R/r
n = 1.0
name = "Poly"
def u0(phi, psi):
    return -9.0
def f(phi, psi, t):
   return 5.0*abs(exp(-0.25*t)*cos(4*phi)*cos(4*psi))
param = param(T, N, M, R, r, tmax, tmin, dt, time, a, name, u0, f)
def my_handler(scene):
    global _frameold
    if (_frameold == scene.frame_current):
        return
    _frameold = scene.frame_current
    param.Recalc(scene.frame_current, False, True)#(фрейм,флаг деформации,флаг изменения цвета)
param.Create()
param.Solution()
param.Recalc(0, False, True)
bpy.app.handlers.frame_change_pre.append(my_handler)