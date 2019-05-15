#ifndef SPACE_H
#define SPACE_H


namespace growth
{

class Space
{
  public:
    Space();

    static const int DEFAULT_DIM;

    static int get_dimension();
    static void set_dimension_(int dim);

  private:
    static int dim_;

  public:
    class Shape
    {
      public:
        Shape();
    };
};
}

#endif // SPACE_H
