#ifndef FILE_OBJECTS
#define FILE_OBJECTS

/* *************************************************************************/
/* copy from netgen  and modify it                                                 */
/* File:   geomobjects.hpp                                                 */
/* Function:  describe fundamental type of geomotry                        */ 
/* Date:   08. jan. 2016                                                   */
/* Class: Poi Vec Mat Bo BoSphere                                      */
/* *************************************************************************/



  template <int D> class Vec;
  template <int D> class Poi;


  template <int D>
  class Poi
  {

  protected:
    double x[D];

  public:
    Poi () { ; }
    Poi (double ax) 
	{
		for (int i = 0; i < D; i++) x[i] = ax; 
	}
    Poi (double ax, double ay) 
    { 
      // static_assert(D==2, "Poi<D> constructor with 2 args called");
      x[0] = ax; x[1] = ay; 
    }
    Poi (double ax, double ay, double az) 
    {
      // static_assert(D==3, "Poi<D> constructor with 3 args called");
      x[0] = ax; x[1] = ay; x[2] = az; 
    }
    Poi (double ax, double ay, double az, double au)
    {
		 // static_assert(D==4, "Poi<D> constructor with 4 args called");
		x[0] = ax; x[1] = ay; x[2] = az; x[3] = au;
	}

    Poi (const Poi<D> & p2)
    {
		for (int i = 0; i < D; i++) x[i] = p2.x[i]; 
	}

    explicit Poi (const Vec<D> & v)
    {
		for (int i = 0; i < D; i++) x[i] = v(i); 
	}


    Poi & operator= (const Poi<D> & p2)
    {
      for (int i = 0; i < D; i++) x[i] = p2.x[i]; 
      return *this;
    }

    Poi & operator= (double val)
    {
      for (int i = 0; i < D; i++) x[i] = val;
      return *this;
    }
	//����ȡС
	const Poi & SetToMin (const Poi<D> & p2)
	{
		for(int i = 0; i < D; i++)
		if (p2.x[i] < x[i]) x[i] = p2.x[i];
		return *this;
	}
	//����ȡ��
	const Poi & SetToMax (const Poi<D> & p2)
	{
		for(int i = 0; i < D; i++)
		if (p2.x[i] > x[i]) x[i] = p2.x[i];
		return *this;
	}

    double & operator() (int i) 
	{ 
		return x[i]; 
	}
    const double & operator() (int i) const
	{
		return x[i]; 
	}
	double & operator[] (int i) 
	{ 
		return x[i]; 
	}
    const double & operator[] (int i) const
	{ 
		return x[i];
	}

    operator const double * () const 
	{ 
		return x;
	}
	//���ĵ�i��ֵ
	void Modify ( int  i , double val )
	{
		x[i] = val; 
	}
	//�������е�ֵ
	void Modify ( double val )
	 { 
		 for (int i = 0; i < D; i++) x[i] = val;
	}

  };


  template <int D>
  class Vec
  {

  protected:
    double x[D];

  public:
    Vec () { ; } // for (int i = 0; i < D; i++) x[i] = 0; }
    Vec (double ax) 
	{ 
		for (int i = 0; i < D; i++) x[i] = ax;
	}
    Vec (double ax, double ay) 
    { 
      // static_assert(D==2, "Vec<D> constructor with 2 args called");
      x[0] = ax; x[1] = ay; 
    }
    Vec (double ax, double ay, double az)
    { 
      // static_assert(D==3, "Vec<D> constructor with 3 args called");
      x[0] = ax; x[1] = ay; x[2] = az; 
    }
    Vec (double ax, double ay, double az, double au)
    { 
		x[0] = ax; x[1] = ay; x[2] = az; x[3] = au; 
	}

    Vec (const Vec<D> & p2)
    { 
		for (int i = 0; i < D; i++) x[i] = p2.x[i]; 
	}

    explicit Vec (const Poi<D> & p)
    { 
		for (int i = 0; i < D; i++) x[i] = p(i); 
	}

    Vec (const Vec<D> & p1, const Vec<D> & p2)
    { 
		for(int i=0; i<D; i++) x[i] = p2(i)-p1(1);
	}
  
    Vec & operator= (const Vec<D> & p2)
    {
      for (int i = 0; i < D; i++) x[i] = p2.x[i]; 
      return *this;
    }

    Vec & operator= (double s)
    {
      for (int i = 0; i < D; i++) x[i] = s;
      return *this;
    }

    double & operator() (int i) 
	{
		return x[i]; 
	}
    const double & operator() (int i) const 
	{ 
		return x[i]; 
	}

    operator const double* () const 
	{ 
		return x; 
	}

	void Modify ( int  i , double val )
	{ 
		x[i] = val; 
	}
	void Modify ( double val )
	{ 
		for (int i = 0; i < D; i++) 
			x[i] = val; 
	}
	//�����ĳ���
    double Length () const
    {
      double l = 0;
      for (int i = 0; i < D; i++)
	      l += x[i] * x[i];
      return sqrt (l);
    }
	//�������ȵ�ƽ��
    double Length2 () const
	{
		double l = 0;
		for (int i = 0; i < D; i++)
			l += x[i] * x[i];
		return l;
	}
	//������һ��
    void Normalize ()
	{
		double l = Length();
		if (l != 0)
			for (int i = 0; i < D; i++)
				x[i] /= l;
	}
	//��ȡ��������
    Vec<D>& GetNormal () const
	{
		Vec<D> vec; 
		double l = Length();
		if (l != 0)
			for (int i = 0; i < D; i++)
				vec.x[i] = x[i] / l;
		return vec;
	}
	
  };




  template <int H, int W=H>
  class Mat
  {
  protected:
    double x[H*W];

  public:
    Mat () { ; }
    Mat (const Mat & b)
    { 
		for (int i = 0; i < H*W; i++) x[i] = b.x[i]; 
	}
  
    Mat & operator= (double s)
    {
      for (int i = 0; i < H*W; i++) x[i] = s;
      return *this;
    }

    Mat & operator= (const Mat & b)
    {
      for (int i = 0; i < H*W; i++) x[i] = b.x[i]; 
      return *this;
    }

    double & operator() (int i, int j) 
	{ 
		return x[i*W+j]; 
	}
    const double & operator() (int i, int j) const 
	{ 
		return x[i*W+j]; 
	}
    double & operator() (int i) 
	{ 
		return x[i]; 
	}
    const double & operator() (int i) const 
	{ 
		return x[i]; 
	}
	//��ȡ��i��
    Vec<H> Col (int i) const
	{
		Vec<H> hv; 
		for (int j = 0; j < H; j++)
			hv(j) = x[j*W+i];
		return hv; 
	}
	//��ȡ��i��
    Vec<W> Row (int i) const
	{
		Vec<W> hv; 
		for (int j = 0; j < W; j++)
			hv(j) = x[i*W+j];
		return hv; 
	}

/*
	����󷽳̻�û��ʵ��
    void Solve (const Vec<H> & rhs, Vec<W> & sol) const
    {
      Mat<W,H> inv;
      CalcInverse (*this, inv);
      sol = inv * rhs;
    }
*/
	void Modify ( int  i , int  j, double val ) 
	{ 
		x[i*W+j] = val; 
	}
	void Modify ( int  i , double val )
	{ 
		x[i] = val; 
	}
	void Modify ( double val )
	{ 
		for (int i = 0; i < H*W; i++) x[i] = val; 
	}
  };

  template <int D>
  class Bo
  {
  protected:
    Poi<D> pmin, pmax;
  public:
    Bo () { ; }

    Bo ( const Poi<D> & p1)
	{
		for (int i = 0; i < D; i++)
			pmin(i) = pmax(i) = p1(i);
	}

    Bo ( const Poi<D> & p1, const Poi<D> & p2)
	{
		for (int i = 0; i < D; i++)
		{
			pmin(i) = min(p1(i), p2(i));
			pmax(i) = max(p1(i), p2(i));
		}
	}

    const Poi<D> & PMin () const 
	{ 
		return pmin; 
	}
    const Poi<D> & PMax () const 
	{ 
		return pmax; 
	}
  
    void Set (const Poi<D> & p)
    { 
		pmin = pmax = p; 
	}

    void Add (const Poi<D> & p)
	{ 
		for (int i = 0; i < D; i++)
		{
			if (p(i) < pmin(i)) pmin(i) = p(i);
			else if (p(i) > pmax(i)) pmax(i) = p(i);
		}
	}

	double & operator[] (int i) 
	{
		if(i<=D) 
			return pmin[i];
		else 
			return pmax[i-D];
	}
	const	double & operator[] (int i) const
	{
		if(i<=D) 
			return pmin[i];
		else 
			return pmax[i-D];
	}
/*
    template <typename T1, typename T2>
    void Set (const IndirectArray<T1, T2> & points)
    {
      Set (points[points.Begin()]);
      for (int i = points.Begin()+1; i < points.End(); i++)
        Add (points[i]);
    }

    template <typename T1, typename T2>
    void Add (const IndirectArray<T1, T2> & points)
    {
      for (int i = points.Begin(); i < points.End(); i++)
        Add (points[i]);
    }
*/
	//���ӵ�����
    Poi<D> Center () const 
	{ 
		Poi<D> c;
		for (int i = 0; i < D; i++)
			c(i) = 0.5 * (pmin(i)+pmax(i)); 
		return c;
	}
	//�Խ��ߵĳ���
    double Diam () const 
	{ 
		return Abs (pmax-pmin); 
	}
	//��ȡ��nr����
    Poi<D> GetPoiNr (int nr) const
	{
		Poi<D> p;
		for (int i = 0; i < D; i++)
		{
			p(i) = (nr & 1) ? pmax(i) : pmin(i);
			nr >>= 1;
		}
		return p;
	}

    bool Intersect (const Bo<D> & box2) const
	{
		for (int i = 0; i < D; i++)
			if (pmin(i) > box2.pmax(i) ||
				pmax(i) < box2.pmin(i)) return 0;
		return 1;
	}
	//�жϵ��Ƿ��ں�����
    bool IsIn (const Poi<D> & p) const
	{
		for (int i = 0; i < D; i++)
			if (p(i) < pmin(i) || p(i) > pmax(i)) return 0;
		return 1;
	}
	//���Ӻ���
    void Increase (double dist)
	{
		for (int i = 0; i < D; i++)
		{
			pmin(i) -= dist;
			pmax(i) += dist;
		}
	}
  };


  template <int D>
  class BoSphere : public Bo<D>
  {
  protected:
	  ///
	  Poi<D> c;
	  ///
	  double diam;
	  ///
	  double inner;
  public:
	  ///
	  BoSphere () { ; }
	  ///
	  BoSphere (const Bo<D> & box) : Bo<D> (box) 
	  { 
		  CalcDiamCenter();
	  }

	  ///
	  BoSphere ( Poi<D> apmin, Poi<D> apmax ): Bo<D> (apmin, apmax)
	  {
		  CalcDiamCenter();
	  }

	  ///
	  const Poi<D> & Center () const 
	  { 
		  return c; 
	  }
	  ///
	  double Diam () const 
	  { 
		  return diam; 
	  }
	  ///
	  double Inner () const 
	  { 
		  return inner; 
	  }


	  ///
	  void GetSubBo (int nr, BoSphere & sbox) const
	  {
		  for (int i = 0; i < D; i++)
		  {
			  if (nr & 1)
			  {
				  sbox.pmin(i) = c(i);
				  sbox.pmax(i) = this->pmax(i);
			  }
			  else
			  {
				  sbox.pmin(i) = this->pmin(i);
				  sbox.pmax(i) = c(i);
			  }
			  sbox.c(i) = 0.5 * (sbox.pmin(i) + sbox.pmax(i));
			  nr >>= 1;
		  }
		  sbox.diam = 0.5 * diam;
		  sbox.inner = 0.5 * inner;
	  }


	  ///
	  void CalcDiamCenter ()
	  {
		  c = Bo<D>::Center ();
		  diam = Dist (this->pmin, this->pmax);

		  inner = this->pmax(0) - this->pmin(0);
		  for (int i = 1; i < D; i++)
			  if (this->pmax(i) - this->pmin(i) < inner)
				  inner = this->pmax(i) - this->pmin(i);
	  }

  };



#endif
