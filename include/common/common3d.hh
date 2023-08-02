#pragma once


template<typename T>
struct vec2 {
    vec2() : x(0), y(0) {}
    vec2(T xIn, T yIn) : x(xIn), y(yIn) {}
    T x;
    T y;
};

template<typename T>
struct vec3 : public vec2<T> {
    vec3() : z(0), vec2<T>(0,0) {}
    vec3(T xIn, T yIn, T zIn) :  z(zIn), vec2<T>(xIn, yIn) {}
    T z;
};

template<typename T>
struct vec4 : public vec3<T> {
    vec4() : vec3<T>(), w(0) {}
    vec4(T xIn, T yIn, T zIn, T wIn) : w(wIn), vec3<T>(xIn, yIn, zIn) {}
    T w;
};

template<typename T>
struct triangle {
    triangle(T i1in, T i2in, T i3in) : i1(i1in), i2(i2in), i3(i3in) {}
    T i1;
    T i2;
    T i3;
};

template<typename T>
struct quad : public triangle<T> {
    quad(T i1in, T i2in, T i3in, T i4in) : triangle<T>(i1in, i2in, i3in), i4(i4in) {}
    T i4;
};
