#ifndef DENSE_TOOLS_H
#define DENSE_TOOLS_H


namespace growth
{

inline int sgn(double x)
{
    return x < 0 ? -1 : (x > 0 ? 1 : 0); 
}

} // namespace growth

#endif /* DENSE_TOOLS_H */
