import numpy as np
from scipy.integrate import odeint

m1 = 43.56
g1 = 9.81
tK = 288
sx = np.pi * 0.077 * 0.077
sy = 0.09637
cx = 0.35
cy = 8.1
zp = 288000 * (1 - (750 / 760) ** (1 / 5.255)) / 6.5
k = 1.4
lambda1 = np.pi / 4
psi = (0 * np.pi) / 6
we = (2 * np.pi) / (24 * 3600)
Vwx = 5
Vwy = 5
alphay = 1.258
alphau = 1.258
v0 = 435
xc = 6000
zc = -100
a = 0


def f(y, t):
    global m1
    global g1
    global tK
    global sx
    global sy
    global cx
    global cy
    global zp
    global k
    global lambda1
    global psi
    global we
    global Vwx
    global Vwy
    global alphay
    global alphau
    j1 = -0.198
    b1 = 1.646
    y1, y2, y3, y4, y5, y6, y7, y8, y9, y10 = y
    return [y2, (-cx * sx * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y5)) * (
                1 - 6.5 * (y5 + zp) / 288000) ** 5.255 * (y2 ** 2 + y6 ** 2 - Vwx * y2) ** (2 + j1 + b1) * y2 / (((
                                                                                                                              k * 8314 * (
                                                                                                                                  tK - 0.006328 * y5) / 28.96) ** 0.5 + Vwx * y2 / np.sqrt(
        y2 ** 2 + y6 ** 2)) ** b1 * (y2 ** 2 + y6 ** 2) ** (0.5 * (3 + j1 + b1))) - (
                        2 * we * (y6 * np.cos(lambda1) * np.cos(psi) - y4 * np.sin(lambda1))), y4,
            (cy * sy * 28.96 * 101325 * np.sign(Vwy)) / (m1 * 8314 * (tK - 0.006328 * y5)) * (
                        1 - 6.5 * (y5 + zp) / 288000) ** 5.255 * (np.abs(Vwy - y4)) ** alphay - (
                        2 * we * (y2 * np.sin(lambda1) - y6 * np.cos(lambda1) * np.sin(psi))), y6,
            -g1 - (cx * sx * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y5)) * (
                        1 - 6.5 * (y5 + zp) / 288000) ** 5.255 * (y2 ** 2 + y6 ** 2 - Vwx * y2) ** (
                        2 + j1 + b1) * y6 / (((k * 8314 * (tK - 0.006328 * y5) / 28.96) ** 0.5 + Vwx * y2 / np.sqrt(
                y2 ** 2 + y6 ** 2)) ** b1 * (y2 ** 2 + y6 ** 2) ** (0.5 * (3 + j1 + b1))) - (
                        2 * we * (y4 * np.cos(lambda1) * np.sin(psi) - y2 * np.cos(lambda1) * np.cos(psi))), y8,
            (np.sign(Vwx) * cy * sy * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y5)) * (
                        1 - 6.5 * (y5 + zp) / 288000) ** 5.255 * (
                np.abs(np.abs(Vwx * y6) / np.sqrt(y2 ** 2 + y6 ** 2) - np.sqrt(y8 ** 2 + y10 ** 2))) ** alphau * np.abs(
                y6) * np.sign(np.abs(Vwx * y6) / np.sqrt(y2 ** 2 + y6 ** 2) - np.sqrt(y8 ** 2 + y10 ** 2)) / np.sqrt(
                y2 ** 2 + y6 ** 2), y10,
            (-np.sign(Vwx) * cy * sy * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y5)) * (
                        1 - 6.5 * (y5 + zp) / 288000) ** 5.255 * (np.abs(
                np.abs(Vwx * y6) / np.sqrt(y2 ** 2 + y6 ** 2) - np.sqrt(y8 ** 2 + y10 ** 2))) ** alphau * y2 * np.sign(
                y6) * (np.sign(np.abs(Vwx * y6) / np.sqrt(y2 ** 2 + y6 ** 2) - np.sqrt(y8 ** 2 + y10 ** 2))) / np.sqrt(
                y2 ** 2 + y6 ** 2)]


y0 = [0, v0 * np.cos(a), 0, 0, 0, v0 * np.sin(a), 0, 0, 0, 0]
t2 = 0
t = np.linspace(0, t2, 500)
w = odeint(f, y0, t)
if zc >= 0:
    while (w[-1][4] + w[-1][8]) <= zc:
        y0 = [0, v0 * np.cos(a), 0, 0, 0, v0 * np.sin(a), 0, 0, 0, 0]
        w = odeint(f, y0, t)
        a = a + 0.1
        while (w[-1][0] + w[-1][6]) <= xc:
            t = np.linspace(0, t2, 500)
            w = odeint(f, y0, t)
            t2 = t2 + 0.1

if zc < 0:
    while (w[-1][0] + w[-1][6]) <= xc:
        t = np.linspace(0, t2, 500)
        w = odeint(f, y0, t)
        t2 = t2 + 0.1
        while (w[-1][4] + w[-1][8]) <= zc:
            y0 = [0, v0 * np.cos(a), 0, 0, 0, v0 * np.sin(a), 0, 0, 0, 0]
            w = odeint(f, y0, t)
            a = a + 0.1

print('KX1 =', w[-1][0])
print('VX1 =', w[-1][1])
print('KY1 =', w[-1][2])
print('VY1 =', w[-1][3])
print('KZ1 =', w[-1][4])
print('VZ1 =', w[-1][5])
print('KU1 =', w[-1][6])
print('VU1 =', w[-1][7])
print('KV1 =', w[-1][8])
print('VV1 =', w[-1][9])
print("Час лету:", t[-1])
print("Кут пострілу:", a)
print('KX1+KU1 =', w[-1][0] + w[-1][6])
print('KZ1+KV1 =', w[-1][4] + w[-1][8])
print("-----------------------------------------")


# Надзвукова швидкість
def h(y, tzv):
    global m1
    global g1
    global tK
    global sx
    global sy
    global cx
    global cy
    global zp
    global k
    global lambda1
    global psi
    global we
    global Vwx
    global Vwy
    global alphay
    global alphau
    j1 = -0.198
    b1 = 1.646
    y1, y2, y3, y4, y5, y6, y7, y8, y9, y10 = y
    return [y2, (-cx * sx * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y5)) * (
                1 - 6.5 * (y5 + zp) / 288000) ** 5.255 * (y2 ** 2 + y6 ** 2 - Vwx * y2) ** (2 + j1 + b1) * y2 / (((
                                                                                                                              k * 8314 * (
                                                                                                                                  tK - 0.006328 * y5) / 28.96) ** 0.5 + Vwx * y2 / np.sqrt(
        y2 ** 2 + y6 ** 2)) ** b1 * (y2 ** 2 + y6 ** 2) ** (0.5 * (3 + j1 + b1))) - (
                        2 * we * (y6 * np.cos(lambda1) * np.cos(psi) - y4 * np.sin(lambda1))), y4,
            (cy * sy * 28.96 * 101325 * np.sign(Vwy)) / (m1 * 8314 * (tK - 0.006328 * y5)) * (
                        1 - 6.5 * (y5 + zp) / 288000) ** 5.255 * (np.abs(Vwy - y4)) ** alphay - (
                        2 * we * (y2 * np.sin(lambda1) - y6 * np.cos(lambda1) * np.sin(psi))), y6,
            -g1 - (cx * sx * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y5)) * (
                        1 - 6.5 * (y5 + zp) / 288000) ** 5.255 * (y2 ** 2 + y6 ** 2 - Vwx * y2) ** (
                        2 + j1 + b1) * y6 / (((k * 8314 * (tK - 0.006328 * y5) / 28.96) ** 0.5 + Vwx * y2 / np.sqrt(
                y2 ** 2 + y6 ** 2)) ** b1 * (y2 ** 2 + y6 ** 2) ** (0.5 * (3 + j1 + b1))) - (
                        2 * we * (y4 * np.cos(lambda1) * np.sin(psi) - y2 * np.cos(lambda1) * np.cos(psi))), y8,
            (np.sign(Vwx) * cy * sy * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y5)) * (
                        1 - 6.5 * (y5 + zp) / 288000) ** 5.255 * (
                np.abs(np.abs(Vwx * y6) / np.sqrt(y2 ** 2 + y6 ** 2) - np.sqrt(y8 ** 2 + y10 ** 2))) ** alphau * np.abs(
                y6) * np.sign(np.abs(Vwx * y6) / np.sqrt(y2 ** 2 + y6 ** 2) - np.sqrt(y8 ** 2 + y10 ** 2)) / np.sqrt(
                y2 ** 2 + y6 ** 2), y10,
            (-np.sign(Vwx) * cy * sy * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y5)) * (
                        1 - 6.5 * (y5 + zp) / 288000) ** 5.255 * (np.abs(
                np.abs(Vwx * y6) / np.sqrt(y2 ** 2 + y6 ** 2) - np.sqrt(y8 ** 2 + y10 ** 2))) ** alphau * y2 * np.sign(
                y6) * (np.sign(np.abs(Vwx * y6) / np.sqrt(y2 ** 2 + y6 ** 2) - np.sqrt(y8 ** 2 + y10 ** 2))) / np.sqrt(
                y2 ** 2 + y6 ** 2)]


# a=0
y0 = [0, v0 * np.cos(a), 0, 0, 0, v0 * np.sin(a), 0, 0, 0, 0]
ts = 0
tzv = np.linspace(0, ts, 500)
ws = odeint(h, y0, tzv)
V1 = np.sqrt(ws[-1][1] ** 2 + ws[-1][5] ** 2)
Vs1 = np.sqrt((k * 8314 * (tK - 0.006328 * ws[-1][4])) / 28.96) + Vwx * ws[-1][1] / np.sqrt(
    ws[-1][1] ** 2 + ws[-1][5] ** 2)
while Vs1 < V1:
    tzv = np.linspace(0, ts, 500)
    ws = odeint(h, y0, tzv)
    V1 = np.sqrt(ws[-1][1] ** 2 + ws[-1][5] ** 2)
    Vs1 = np.sqrt((k * 8314 * (tK - 0.006328 * ws[-1][4])) / 28.96) + Vwx * ws[-1][1] / np.sqrt(
        ws[-1][1] ** 2 + ws[-1][5] ** 2)
    ts = ts + 0.1
print("Vs1:", Vs1)
print("V1:", V1)
print("Дозвуковий час tzv:", tzv[-1])
print("KX2:", ws[-1][0])
print('VX2 =', ws[-1][1])
print('KY2 =', ws[-1][2])
print('VY2 =', ws[-1][3])
print('KZ2 =', ws[-1][4])
print('VZ2 =', ws[-1][5])
print('KU2 =', ws[-1][6])
print('VU2 =', ws[-1][7])
print('KV2 =', ws[-1][8])
print('VV2 =', ws[-1][9])

# Другий етап
print("-----------------------------------------")


def g(y, t100):
    global m1
    global g1
    global tK
    global sx
    global sy
    global cx
    global cy
    global zp
    global k
    global lambda1
    global psi
    global we
    global Vwx
    global Vwy
    global alphay
    global alphau
    j2 = -0.09
    b2 = 6.548
    y10, y20, y30, y40, y50, y60, y70, y80, y90, y100 = y
    return [y20, (-cx * sx * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y50)) * (
                1 - 6.5 * (y50 + zp) / 288000) ** 5.255 * (y20 ** 2 + y60 ** 2 - Vwx * y20) ** (2 + j2 + b2) * y20 / (((
                                                                                                                                   k * 8314 * (
                                                                                                                                       tK - 0.006328 * y50) / 28.96) ** 0.5 + Vwx * y20 / np.sqrt(
        y20 ** 2 + y60 ** 2)) ** b2 * (y20 ** 2 + y60 ** 2) ** (0.5 * (3 + j2 + b2))) - (
                        2 * we * (y60 * np.cos(lambda1) * np.cos(psi) - y40 * np.sin(lambda1))), y40,
            (cy * sy * 28.96 * 101325 * np.sign(Vwy)) / (m1 * 8314 * (tK - 0.006328 * y50)) * (
                        1 - 6.5 * (y50 + zp) / 288000) ** 5.255 * (np.abs(Vwy - y40)) ** alphay - (
                        2 * we * (y20 * np.sin(lambda1) - y60 * np.cos(lambda1) * np.sin(psi))), y60,
            -g1 - (cx * sx * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y50)) * (
                        1 - 6.5 * (y50 + zp) / 288000) ** 5.255 * (y20 ** 2 + y60 ** 2 - Vwx * y20) ** (
                        2 + j2 + b2) * y60 / (((k * 8314 * (tK - 0.006328 * y50) / 28.96) ** 0.5 + Vwx * y20 / np.sqrt(
                y20 ** 2 + y60 ** 2)) ** b2 * (y20 ** 2 + y60 ** 2) ** (0.5 * (3 + j2 + b2))) - (
                        2 * we * (y40 * np.cos(lambda1) * np.sin(psi) - y20 * np.cos(lambda1) * np.cos(psi))), y80,
            (np.sign(Vwx) * cy * sy * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y50)) * (
                        1 - 6.5 * (y50 + zp) / 288000) ** 5.255 * (np.abs(
                np.abs(Vwx * y60) / np.sqrt(y20 ** 2 + y60 ** 2) - np.sqrt(y80 ** 2 + y100 ** 2))) ** alphau * np.abs(
                y60) * np.sign(
                np.abs(Vwx * y60) / np.sqrt(y20 ** 2 + y60 ** 2) - np.sqrt(y80 ** 2 + y100 ** 2)) / np.sqrt(
                y20 ** 2 + y60 ** 2), y100,
            (-np.sign(Vwx) * cy * sy * 28.96 * 101325) / (m1 * 8314 * (tK - 0.006328 * y50)) * (
                        1 - 6.5 * (y50 + zp) / 288000) ** 5.255 * (np.abs(
                np.abs(Vwx * y60) / np.sqrt(y20 ** 2 + y60 ** 2) - np.sqrt(
                    y80 ** 2 + y100 ** 2))) ** alphau * y20 * np.sign(y60) * (
                np.sign(np.abs(Vwx * y60) / np.sqrt(y20 ** 2 + y60 ** 2) - np.sqrt(y80 ** 2 + y100 ** 2))) / np.sqrt(
                y20 ** 2 + y60 ** 2)]


y01 = [ws[-1][0], ws[-1][1], ws[-1][2], ws[-1][3], ws[-1][4], ws[-1][5], ws[-1][6], ws[-1][7], ws[-1][8], ws[-1][9]]
t3 = 0
t100 = np.linspace(0, t3, 500)
w1 = odeint(g, y01, t100)
# a=0
"""if zc>=0:     
    while (w1[-1][4] + w1[-1][8]) <= zc:
        y01 = [ws[-1][0], ws[-1][1], ws[-1][2], ws[-1][3], ws[-1][4], ws[-1][5], ws[-1][6], ws[-1][7], ws[-1][8], ws[-1][9]]
        w1 = odeint(g, y01, t100)
        a = a + 0.1
        while (w1[-1][0] + w1[-1][6]) <= xc:
            t100 = np.linspace(0, t3, 500)
            w1 = odeint(g, y01, t100)
            t3 = t3 + 0.1 """

if zc > 0:
    while (w1[-1][0] + w1[-1][6]) <= xc:
        t100 = np.linspace(0, t3, 500)
        w1 = odeint(g, y01, t100)
        t3 = t3 + 0.1
        while (w1[-1][4] + w1[-1][8]) <= zc:
            y01 = [ws[-1][0], ws[-1][1], ws[-1][2], ws[-1][3], ws[-1][4], ws[-1][5], ws[-1][6], ws[-1][7], ws[-1][8],
                   ws[-1][9]]
            w1 = odeint(g, y01, t100)
            a = a + 0.1
        print("Кут пострілу:", a)
print('KX2 =', w1[-1][0])
print('VX2 =', w1[-1][1])
print('KY2 =', w1[-1][2])
print('VY2 =', w1[-1][3])
print('KZ2 =', w1[-1][4])
print('VZ2 =', w1[-1][5])
print('KU2 =', w1[-1][6])
print('VU2 =', w1[-1][7])
print('KV2 =', w1[-1][8])
print('VV2 =', w1[-1][9])
print("Кінцевий час лету:", t100[-1] + tzv[-1])
print("Кут пострілу:", a)





































