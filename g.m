function [g_ijk] = g(rij,rik,rjk,c,d,h)

cosTh = (rij^2 + rik^2 - rjk^2)/(2 * rij * rik);

g_ijk = 1+ c^2/d^2 - c^2/(d^2+(h-cosTh)^2);