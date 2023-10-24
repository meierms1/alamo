// 
// Initialize a field using a bitmap image. (2D only)
//
// Note that in GIMP, you must select "do not write color space information"
// and "24 bit R8 G8 B8" when exporting the BMP file.
//

#ifndef IC_EulerAngles_H
#define IC_EulerAngles_H
#include <cmath>

#include "IC/IC.H"
#include "Util/Util.H"
#include "Util/BMP.H"
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace IC
{
class EulerAngles : public IC
{
public:
    //enum Type {XYZ, XY, YZ, XZ};
    enum Fit {Stretch,FitWidth,FitHeight};
    enum Channel {R=0, G=1, B=2};

    EulerAngles (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom)  {}

    EulerAngles (amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name) : EulerAngles(_geom)
    {
        pp.queryclass(name,*this);
    }

    void Define(std::string eulerfilename)
    {
        bmp.Define(eulerfilename);//"Interface_Blur2.bmp");
    }
  
    void Add(const int &lev, Set::Field<Set::Scalar> &a_field, Set::Scalar)
    {        
        Util::Assert(INFO,TEST(a_field[0]->nComp() == 3));
        Set::Vector DX(geom[lev].CellSize());
        amrex::Box domain = geom[lev].Domain();

        //Set::Scalar width  = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];
        //Set::Scalar height = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];

        Set::Scalar img_width = (Set::Scalar)(bmp.nx-1);
        Set::Scalar img_height = (Set::Scalar)(bmp.ny-1);
        Set::Scalar img_dx = 1.0;
        Set::Scalar img_dy = 1.0;

        amrex::IndexType type = a_field[lev]->ixType();
        
        domain.convert(type);

        Set::Vector domlo(AMREX_D_DECL(geom[lev].ProbLo()[0],geom[lev].ProbLo()[1],0.0));
        Set::Vector domhi(AMREX_D_DECL(geom[lev].ProbHi()[0],geom[lev].ProbHi()[1],0.0));


        for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();
            bx = bx & domain;


            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            //amrex::Array4<Set::Scalar> const& fieldR = a_fieldR[lev]->array(mfi);   //newly added line
            //amrex::Array4<Set::Scalar> const& fieldG = a_fieldG[lev]->array(mfi);  //newly added line
            //amrex::Array4<Set::Scalar> const& fieldB = a_fieldB[lev]->array(mfi); //newly added line
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                Set::Vector x = Set::Vector::Zero();
                // NODE
                if (type == amrex::IndexType::TheNodeType())
                {
                    x(0) = domlo(0) + ((amrex::Real)(i)) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + ((amrex::Real)(j)) * geom[lev].CellSize()[1];
                }
                else if (type == amrex::IndexType::TheCellType())
                {
                    x(0) = domlo(0) + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
                }

                //Set::Scalar x, y;
                Set::Vector ximg;

                if (fit == Fit::Stretch)
                {
                    ximg(0) = (x(0) - domlo(0)) / (domhi(0) - domlo(0)); 
                    ximg(1) = (x(1) - domlo(1)) / (domhi(1) - domlo(1)); 
                }
                else if (fit == Fit::FitWidth)
                {
                    Set::Scalar aspect_ratio = img_width / img_height;
                    ximg(0) = (x(0) - domlo(0)) / (domhi(0) - domlo(0)); 
                    ximg(1) = (x(1) - domlo(1)) / (domhi(1) - domlo(1)); 
                    ximg(1) -= 0.5 - 0.5 / aspect_ratio;
                    ximg(1) *= aspect_ratio;
                }
                else if (fit == Fit::FitHeight)
                {
                    Set::Scalar aspect_ratio = img_height / img_width;
                    ximg(0) = (x(0) - domlo(0)) / (domhi(0) - domlo(0)); 
                    ximg(1) = (x(1) - domlo(1)) / (domhi(1) - domlo(1)); 
                    ximg(0) -= 0.5 - 0.5 / aspect_ratio;
                    ximg(0) *= aspect_ratio;
                }

                ximg(0) = std::min(ximg(0),1.0); ximg(1) = std::min(ximg(1),1.0);
                ximg(0) = std::max(ximg(0),0.0); ximg(1) = std::max(ximg(1),0.0);
                
                ximg(0) *= img_width;
                ximg(1) *= img_height;

                int I = (int)(ximg(0));
                int J = (int)(ximg(1));

                Set::Scalar x1 = I*img_dx, x2 = (I+1)*img_dx;
                Set::Scalar y1 = J*img_dy, y2 = (J+1)*img_dy;

                if (I < bmp.nx-1 && J < bmp.ny - 1)
                {
                    //Set::Scalar fQ11R = ((Set::Scalar)(bmp(I,  J  )[channel]) - min) / (max - min);
                    //Set::Scalar fQ12R = ((Set::Scalar)(bmp(I,  J+1)[channel]) - min) / (max - min);
                    //Set::Scalar fQ21R = ((Set::Scalar)(bmp(I+1,J  )[channel]) - min) / (max - min);
                    //Set::Scalar fQ22R = ((Set::Scalar)(bmp(I+1,J+1)[channel]) - min) / (max - min);
                    Set::Scalar fQ11R = ((Set::Scalar)(bmp(I,  J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ12R = ((Set::Scalar)(bmp(I,  J+1)[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ21R = ((Set::Scalar)(bmp(I+1,J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ22R = ((Set::Scalar)(bmp(I+1,J+1)[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ11G = ((Set::Scalar)(bmp(I,  J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ12G = ((Set::Scalar)(bmp(I,  J+1)[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ21G = ((Set::Scalar)(bmp(I+1,J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ22G = ((Set::Scalar)(bmp(I+1,J+1)[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ11B = ((Set::Scalar)(bmp(I,  J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ12B = ((Set::Scalar)(bmp(I,  J+1)[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ21B = ((Set::Scalar)(bmp(I+1,J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ22B = ((Set::Scalar)(bmp(I+1,J+1)[channel3]) - min[2]) / (max[2] - min[2]);
                    

                    //field(i,j,k) =
                    //    (fQ11*(x2-ximg(0))*(y2-ximg(1)) + fQ21*(ximg(0)-x1)*(y2-ximg(1)) + fQ12*(x2-ximg(0))*(ximg(1)-y1) + fQ22*(ximg(0)-x1)*(ximg(1)-y1)) / (img_dx * img_dy);
                    
                    field(i,j,k,0) =
                        (fQ11R*(x2-ximg(0))*(y2-ximg(1)) + fQ21R*(ximg(0)-x1)*(y2-ximg(1)) + fQ12R*(x2-ximg(0))*(ximg(1)-y1) + fQ22R*(ximg(0)-x1)*(ximg(1)-y1)) / (img_dx * img_dy);
                        
                    field(i,j,k,1) =
                        (fQ11G*(x2-ximg(0))*(y2-ximg(1)) + fQ21G*(ximg(0)-x1)*(y2-ximg(1)) + fQ12G*(x2-ximg(0))*(ximg(1)-y1) + fQ22G*(ximg(0)-x1)*(ximg(1)-y1)) / (img_dx * img_dy);
                    
                    field(i,j,k,2) =
                        (fQ11B*(x2-ximg(0))*(y2-ximg(1)) + fQ21B*(ximg(0)-x1)*(y2-ximg(1)) + fQ12B*(x2-ximg(0))*(ximg(1)-y1) + fQ22B*(ximg(0)-x1)*(ximg(1)-y1)) / (img_dx * img_dy);
                    
                    
                    
               
                }
                else if (I == bmp.nx-1 && J <  bmp.ny - 1)
                {
                    //Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[channel]) - min) / (max - min);
                    //Set::Scalar fQ12 = ((Set::Scalar)(bmp(I  ,J+1)[channel]) - min) / (max - min);
                    //field(i,j,k) = fQ11 + (fQ12-fQ11) * (ximg(1) - y1);
                    
                    Set::Scalar fQ11R = ((Set::Scalar)(bmp(I  ,J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ12R = ((Set::Scalar)(bmp(I  ,J+1)[channel1]) - min[0]) / (max[0] - min[0]);
                    field(i,j,k,0) = fQ11R + (fQ12R-fQ11R) * (ximg(1) - y1);
                    
                    Set::Scalar fQ11G = ((Set::Scalar)(bmp(I  ,J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ12G = ((Set::Scalar)(bmp(I  ,J+1)[channel2]) - min[1]) / (max[1] - min[1]);
                    field(i,j,k,1) = fQ11G + (fQ12G-fQ11G) * (ximg(1) - y1);
                    
                    Set::Scalar fQ11B = ((Set::Scalar)(bmp(I  ,J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ12B = ((Set::Scalar)(bmp(I  ,J+1)[channel3]) - min[2]) / (max[2] - min[2]);
                    field(i,j,k,2) = fQ11B + (fQ12B-fQ11B) * (ximg(1) - y1);
                }
                else if (I <  bmp.nx-1 && J == bmp.ny - 1)
                {
                    //Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[channel]) - min) / (max - min);
                    //Set::Scalar fQ21 = ((Set::Scalar)(bmp(I+1,J  )[channel]) - min) / (max - min);
                    //field(i,j,k) = fQ11 + (fQ21-fQ11) * (ximg(0) - x1);
                    
                    Set::Scalar fQ11R = ((Set::Scalar)(bmp(I  ,J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    Set::Scalar fQ21R = ((Set::Scalar)(bmp(I+1,J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    field(i,j,k,0) = fQ11R + (fQ21R-fQ11R) * (ximg(0) - x1);
                    
                    Set::Scalar fQ11G = ((Set::Scalar)(bmp(I  ,J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    Set::Scalar fQ21G = ((Set::Scalar)(bmp(I+1,J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    field(i,j,k,1) = fQ11G + (fQ21G-fQ11G) * (ximg(0) - x1);
                    
                    Set::Scalar fQ11B = ((Set::Scalar)(bmp(I  ,J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    Set::Scalar fQ21B = ((Set::Scalar)(bmp(I+1,J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    field(i,j,k,2) = fQ11B + (fQ21B-fQ11B) * (ximg(0) - x1);
                    
                }
                else if (I ==  bmp.nx-1 && J == bmp.ny - 1)
                {
                    //Set::Scalar fQ11 = ((Set::Scalar)(bmp(I  ,J  )[channel]) - min) / (max - min);
                    //field(i,j,k) = fQ11;
                    Set::Scalar fQ11R = ((Set::Scalar)(bmp(I  ,J  )[channel1]) - min[0]) / (max[0] - min[0]);
                    field(i,j,k,0) = fQ11R;
                    Set::Scalar fQ11G = ((Set::Scalar)(bmp(I  ,J  )[channel2]) - min[1]) / (max[1] - min[1]);
                    field(i,j,k,1) = fQ11G;
                    Set::Scalar fQ11B = ((Set::Scalar)(bmp(I  ,J  )[channel3]) - min[2]) / (max[2] - min[2]);
                    field(i,j,k,2) = fQ11B;
                }
                else
                {
                    //field(i,j,k) = 0.0;
                    field(i,j,k,0) = 0.0;
                    field(i,j,k,1) = 0.0;
                    field(i,j,k,2) = 0.0;
                }

                if (field.nComp() > 1) field(i,j,k,1) = 1.0 - field(i,j,k,0);
                //if (field.nComp() > 1) field(i,j,k,1) = 1.0 - field(i,j,k,0);
                //if (field.nComp() > 1) field(i,j,k,1) = 1.0 - field(i,j,k,0);
                //if (field.nComp() > 1) field(i,j,k,1) = 1.0 - field(i,j,k,0);

            });
        }
        a_field[lev]->FillBoundary();
        //a_field[lev]->FillBoundary();
        //a_fieldG[lev]->FillBoundary();
        //a_fieldB[lev]->FillBoundary();
    };
  
private:
    //Util::BMP bmp;
    Util::BMP bmp;
    Fit fit = Fit::Stretch;
    Channel channel1 = Channel::R;
    Channel channel2 = Channel::G;
    Channel channel3 = Channel::B;
    //Set::Scalar min=0.0, max=255.0;
    Set::Scalar min[3] = {0.0, 0.0, 0.0};
    Set::Scalar max[3] = {255.0, 255.0, 255.0};

public:
    static void Parse(EulerAngles & value, IO::ParmParse & pp)    //Replaced BMP with EulerAngles
    {
        std::string filename;
        pp.query("filename",filename); // BMP filename.
        value.bmp.Define(filename);

        std::string fit = "stretch";
        pp.query("fit",fit); // How to fit. (options: stretch, fitheight, fitwidth)
        if (fit=="stretch")        value.fit = Fit::Stretch;
        else if (fit=="fitheight") value.fit = Fit::FitHeight;
        else if (fit=="fitwidth")  value.fit = Fit::FitWidth;
        else Util::Abort(INFO,"Invalid value for bmp fit - should be stretch/fitheight/fitwidth but received '",fit,"'");

        //std::string channel = "g";
        //pp.query("channel",channel); // Color channel to use (options: r, R, g, G, b, B)
        //if      (channel=="r" || channel=="R") value.channel = Channel::R;
        //else if (channel=="g" || channel=="G") value.channel = Channel::G;
        //else if (channel=="b" || channel=="B") value.channel = Channel::B;
        //else Util::Abort(INFO,"Invalid value for bmp channel - should be r/g/b but received '",channel,"'");
        //pp.query("channel1",value.channel1);
        //pp.query("channel2",value.channel2);
        //pp.query("channel3",value.channel3);
        value.channel1 = Channel::R;
        value.channel2 = Channel::G;
        value.channel3 = Channel::B;

        //value.min = (Set::Scalar) value.bmp.min()[value.channel];
        //value.max = (Set::Scalar) value.bmp.max()[value.channel];
        //value.min() = (Set::Scalar) value.bmp.min()[value.channel1];
        //value.max() = (Set::Scalar) value.bmp.max()[value.channel1];
        //value.min() = (Set::Scalar) value.bmp.min()[value.channel2];
        //value.max() = (Set::Scalar) value.bmp.max()[value.channel2];
        //value.minB( = (Set::Scalar) value.bmp.min()[value.channel3];
        //value.maxB = (Set::Scalar) value.bmp.max()[value.channel3];
        
        min[0] = (Set::Scalar) value.bmp.min()[value.channel1];
        max[0] = (Set::Scalar) value.bmp.max()[value.channel1];

        min[1] = (Set::Scalar) value.bmp.min()[value.channel2];
        max[1] = (Set::Scalar) value.bmp.max()[value.channel2];

        min[2] = (Set::Scalar) value.bmp.min()[value.channel3];
        max[2] = (Set::Scalar) value.bmp.max()[value.channel3];

        
        //pp.query("min",value.min); // Scaling value - minimum (default: 0.0)
        //pp.query("max",value.max); // Scaling value - maximum (default: 255.0)
        pp.query("minR",min[0]);
        pp.query("maxR",max[0]);
        pp.query("minG",min[1]);
        pp.query("maxG",max[1]);
        pp.query("minB",min[2]);
        pp.query("maxB",max[2]);
}    
};
}
#endif
