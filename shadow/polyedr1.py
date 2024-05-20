from math import pi, sqrt
from functools import reduce
from operator import add
from common.r3 import R3
from common.tk_drawer import TkDrawer


def square(n, k):
    x0 = k[0][0]
    y0 = k[0][1]
    x1 = k[1][0]
    y1 = k[1][1]
    sum1 = x0 * y1
    sum2 = y0 * x1
    for i in range(2, n):
        sum1 = sum1 + x1 * k[i][1]
        sum2 = sum2 + y1 * k[i][0]
        x1 = k[i][0]
        y1 = k[i][1]
    sum1 = sum1 + x1 * y0
    sum2 = sum2 + x0 * y1
    return abs((sum1 - sum2) / 2)


class Segment:
    """ Одномерный отрезок """
    # Параметры конструктора: начало и конец отрезка (числа)

    def __init__(self, beg, fin):
        self.beg, self.fin = beg, fin

    # Отрезок вырожден?
    def is_degenerate(self):
        return self.beg >= self.fin

    # Пересечение с отрезком
    def intersect(self, other):
        if other.beg > self.beg:
            self.beg = other.beg
        if other.fin < self.fin:
            self.fin = other.fin
        return self

    # Разность отрезков
    # Разность двух отрезков всегда является списком из двух отрезков!
    def subtraction(self, other):
        return [Segment(
            self.beg, self.fin if self.fin < other.beg else other.beg),
            Segment(self.beg if self.beg > other.fin else other.fin, self.fin)]


class Edge:
    """ Ребро полиэдра """
    # Начало и конец стандартного одномерного отрезка
    SBEG, SFIN = 0.0, 1.0

    # Параметры конструктора: начало и конец ребра (точки в R3)
    def __init__(self, beg, fin):
        self.beg, self.fin = beg, fin
        # Список «просветов»
        self.gaps = [Segment(Edge.SBEG, Edge.SFIN)]

    # Учёт тени от одной грани
    def shadow(self, facet):
        # «Вертикальная» грань не затеняет ничего
        if facet.is_vertical():
            return
        # Нахождение одномерной тени на ребре
        shade = Segment(Edge.SBEG, Edge.SFIN)
        for u, v in zip(facet.vertexes, facet.v_normals()):
            shade.intersect(self.intersect_edge_with_normal(u, v))
            if shade.is_degenerate():
                return

        shade.intersect(
            self.intersect_edge_with_normal(
                facet.vertexes[0], facet.h_normal()))
        if shade.is_degenerate():
            return
        # Преобразование списка «просветов», если тень невырождена
        gaps = [s.subtraction(shade) for s in self.gaps]
        self.gaps = [
            s for s in reduce(add, gaps, []) if not s.is_degenerate()]

    # Преобразование одномерных координат в трёхмерные
    def r3(self, t):
        return self.beg * (Edge.SFIN - t) + self.fin * t

    # Пересечение ребра с полупространством, задаваемым точкой (a)
    # на плоскости и вектором внешней нормали (n) к ней
    def intersect_edge_with_normal(self, a, n):
        f0, f1 = n.dot(self.beg - a), n.dot(self.fin - a)
        if f0 >= 0.0 and f1 >= 0.0:
            return Segment(Edge.SFIN, Edge.SBEG)
        if f0 < 0.0 and f1 < 0.0:
            return Segment(Edge.SBEG, Edge.SFIN)
        x = - f0 / (f1 - f0)
        return Segment(Edge.SBEG, x) if f0 < 0.0 else Segment(x, Edge.SFIN)


class Facet:
    """ Грань полиэдра """
    # Параметры конструктора: список вершин

    def __init__(self, vertexes):
        self.vertexes = vertexes

    # «Вертикальна» ли грань?
    def is_vertical(self):
        return self.h_normal().dot(Polyedr.V) == 0.0

    # Нормаль к «горизонтальному» полупространству
    def h_normal(self):
        n = (
            self.vertexes[1] - self.vertexes[0]).cross(
            self.vertexes[2] - self.vertexes[0])
        return n * (-1.0) if n.dot(Polyedr.V) < 0.0 else n

    # Нормали к «вертикальным» полупространствам, причём k-я из них
    # является нормалью к грани, которая содержит ребро, соединяющее
    # вершины с индексами k-1 и k
    def v_normals(self):
        return [self._vert(x) for x in range(len(self.vertexes))]

    # Вспомогательный метод
    def _vert(self, k):
        n = (self.vertexes[k] - self.vertexes[k - 1]).cross(Polyedr.V)
        return n * \
            (-1.0) if n.dot(self.vertexes[k - 1] - self.center()) < 0.0 else n

    # Центр грани
    def center(self):
        return sum(self.vertexes, R3(0.0, 0.0, 0.0)) * \
            (1.0 / len(self.vertexes))


class Polyedr:
    """ Полиэдр """
    # вектор проектирования
    V = R3(0.0, 0.0, 1.0)

    # Параметры конструктора: файл, задающий полиэдр
    def __init__(self, file):

        # списки вершин, рёбер и граней полиэдра
        self.vertexes, self.edges, self.facets, self.summa = [], [], [], 0

        # список строк файла
        with (open(file) as f):
            for i, line in enumerate(f):
                if i == 0:
                    # обрабатываем первую строку; buf - вспомогательный массив
                    buf = line.split()
                    # коэффициент гомотетии
                    c = float(buf.pop(0))
                    # углы Эйлера, определяющие вращение
                    alpha, beta, gamma = (float(x) * pi / 180.0 for x in buf)
                elif i == 1:
                    # во второй строке число вершин, граней и рёбер полиэдра
                    nv, nf, ne = (int(x) for x in line.split())
                elif i < nv + 2:
                    # задание всех вершин полиэдра
                    x, y, z = (float(x) for x in line.split())
                    if abs(x) >= 1 or \
                        abs(y) >= 1 or \
                            abs(x) <= 0.5 and abs(y) <= 0.5:
                        self.vertexes.append((R3(x, y, z).rz(
                            alpha).ry(beta).rz(gamma) * c, False, R3(x, y, z)))
                    else:
                        self.vertexes.append((R3(x, y, z).rz(
                            alpha).ry(beta).rz(gamma) * c, True, R3(x, y, z)))
                else:
                    mas1 = []
                    mas2 = []
                    mas3 = []
                    needle = True
                    sumx = 0
                    sumy = 0
                    # вспомогательный массив
                    buf = line.split()
                    # количество вершин очередной грани
                    size = int(buf.pop(0))
                    # массив вершин этой грани
                    for n in buf:
                        if not self.vertexes[int(n) - 1][1]:
                            needle = False
                        vertexes = [self.vertexes[int(n) - 1][0] for n in buf]
                        vertexes2 = [self.vertexes[int(n) - 1][2] for n in buf]
                        sumx += self.vertexes[int(n) - 1][2].x
                        sumy += self.vertexes[int(n) - 1][2].y
                    if needle:
                        if abs(sumx / size) >= 1 or \
                            abs(sumy / size) >= 1 or \
                                abs(sumx / size) <= 0.5 and \
                                abs(sumy / size) <= 0.5:
                            needle = False
                    if needle:
                        for i in range(size):
                            element1 = (vertexes2[i].x, vertexes2[i].y)
                            mas1.append(element1)
                            element2 = (vertexes2[i].x, vertexes2[i].z)
                            mas2.append(element2)
                            element3 = (vertexes2[i].y, vertexes2[i].z)
                            mas3.append(element3)
                        xv = vertexes2[1].y * vertexes2[2].z + \
                            vertexes2[2].y * vertexes2[3].z + \
                            vertexes2[3].y * vertexes2[1].z - \
                            vertexes2[1].z * vertexes2[2].y - \
                            vertexes2[2].z * vertexes2[3].y - \
                            vertexes2[3].z * vertexes2[1].y
                        yv = vertexes2[1].z * vertexes2[2].x + \
                            vertexes2[2].z * vertexes2[3].x + \
                            vertexes2[3].z * vertexes2[1].x - \
                            vertexes2[1].x * vertexes2[2].z - \
                            vertexes2[2].x * vertexes2[3].z - \
                            vertexes2[3].x * vertexes2[1].z
                        zv = vertexes2[1].x * vertexes2[2].y + \
                            vertexes2[2].x * vertexes2[3].y + \
                            vertexes2[3].x * vertexes2[1].y - \
                            vertexes2[1].y * vertexes2[2].x - \
                            vertexes2[2].y * vertexes2[3].x - \
                            vertexes2[3].y * vertexes2[1].x
                        angle1 = abs(zv / sqrt(xv ** 2 + yv ** 2 + zv ** 2))
                        angle2 = abs(xv / sqrt(xv ** 2 + yv ** 2 + zv ** 2))
                        angle3 = abs(yv / sqrt(xv ** 2 + yv ** 2 + zv ** 2))
                        ans1 = square(size, mas1)
                        ans2 = square(size, mas3)
                        ans3 = square(size, mas2)
                        if angle1 != 0:
                            self.summa = self.summa + ans1 / angle1
                        elif angle2 != 0:
                            self.summa = self.summa + ans2 / angle2
                        else:
                            self.summa = self.summa + ans3 / angle3
                    # задание рёбер грани
                    for n in range(size):
                        self.edges.append(Edge(vertexes[n - 1], vertexes[n]))
                    # задание самой грани
                    self.facets.append(Facet(vertexes))

    # Метод изображения полиэдра
    def draw(self, tk):  # pragma: no cover
        tk.clean()
        for e in self.edges:
            for f in self.facets:
                e.shadow(f)
            for s in e.gaps:
                tk.draw_line(e.r3(s.beg), e.r3(s.fin))
        print('Площадь граней удолетворяющих условию = ', self.summa)
