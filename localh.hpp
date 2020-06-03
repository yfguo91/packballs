#ifndef LOCALH
#define LOCALH

/**************************************************************************/
/* File:   localh.hh                                                      */
/* Function: Control of 3D mesh grading and restrict sphere radius        */
/* Date:   24. jan. 2016                                                   */
/**************************************************************************/

#include<vector>
#include<ostream>
#include"geomobjects.hpp"
using namespace std;

  /// box for grading
  class GradingBo
  {
    /// xmid
    double xmid[3];
    /// half edgelength
    double h2;
    ///
    GradingBo * childs[8];
    ///
    GradingBo * father;
    ///
    double hopt;
    ///
  public:

    //struct 
    //{
    //  unsigned int cutboundary:1;
    //  unsigned int isinner:1;
    //  unsigned int oldcell:1;
    //  unsigned int pinner:1;
    //} flags;

    ///
    GradingBo (const double * ax1, const double * ax2);
    ///
    void DeleteChilds();
    ///

    Poi<3> PMid() const 
	{ 
		return Poi<3> (xmid[0], xmid[1], xmid[2]); 
	}
    double H2() const 
	{ 
		return h2; 
	}

    friend class LocalH;

  };




  /**
     Control of 3D mesh grading
  */
  class LocalH 
  {
    ///
    GradingBo * root;
    ///
    double grading;
    ///
    vector<GradingBo*> boxes;
    ///
    Bo<3> boundingbox;
  public:
    ///
    LocalH (const Poi<3> & pmin, const Poi<3> & pmax, double grading);
    ///
    LocalH (const Bo<3> & box, double grading);
    ///
    ~LocalH();
    ///
    void Delete();
    ///
    void SetGrading (double agrading)
	{ 
		grading = agrading;
	}
    ///
    void SetH (const Poi<3> & x, double h);
    ///
    double GetH (const Poi<3> & x) const;
    /// minimal h in box (pmin, pmax)
    double GetMinH (const Poi<3> & pmin, const Poi<3> & pmax) const;

    /// mark boxes intersecting with boundary-box
    // void CutBoundary (const Poi3d & pmin, const Poi3d & pmax)
    // { CutBoundaryRec (pmin, pmax, root); }
    void CutBoundary (const Bo<3> & box)
    { 
		CutBoundaryRec (box.PMin(), box.PMax(), root); 
	}
  

    /// widen refinement zone
    void WidenRefinement ();

    int GetNBoes () 
	{
		return boxes.size();
	} 
    const Bo<3> & GetBoundingBo () const
    { 
		return boundingbox; 
	}
    ///
    void PrintMemInfo (ostream & ost) const;
  private:
    /// 
    double GetMinHRec (const Poi<3> & pmin, const Poi<3> & pmax,const GradingBo * box) const;
    ///
    void CutBoundaryRec (const Poi<3> & pmin, const Poi<3> & pmax,GradingBo * box);

    friend ostream & operator<< (ostream & ost, const LocalH & loch);
  };




  inline ostream & operator<< (ostream & ost, const GradingBo & box)
  {
    ost << "gradbox, pmid = " << box.PMid() << ", h2 = " << box.H2() 
	<< endl;
    return ost;
  }

  inline ostream & operator<< (ostream & ost, const LocalH & loch)
  {
    for (int i = 0; i < loch.boxes.size(); i++)
      ost << "box[" << i << "] = " << *(loch.boxes[i]);
    return ost;
  }


#endif
