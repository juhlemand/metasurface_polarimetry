#include "colors.inc"
#include "shapes.inc"
#include "textures.inc"

#declare height_array = array[5] {1.0, 3.0, 2.0, 1.0, 2.0};
#declare width_array = array[5] {2.0, 3.0, 1.0, 3.5, 2.3};
#declare separation = 4.0;
#declare thick = 2;

#declare num_rows = 20;
#declare num_columns = 4;
#declare column_width = 5*separation;


#declare laser_filter = 0.25; // defines transparency of laser diffraction orders
#declare move_forward  = 0; // how far arrows and polarization ellipses are translated in z
#declare ellipse_displacement = 1.45; // dictates how far pol ellipses lie from arrow tips
#declare sphere_rad = 0.4; // dictates how thick the pol ellipses are
#declare ellipse_rad = 8; // base size of pol ellipse
#declare arrow_height = 0.55*ellipse_rad; // height of arrows in pol ellipses
#declare arrow_rad = 0.15*ellipse_rad; // radius of base of arrows in pol ellipses

background { White }
#declare row = 0;
#while (row < num_rows)
    #declare column = 0;
    #while (column < num_columns)
        #declare Index1 = 0;
        #while(Index1<5)
          
          
            box {
                <separation*(Index1)-width_array[Index1]/2 + column*column_width, -height_array[Index1]/2+separation*row, thick>, <separation*(Index1)+width_array[Index1]/2 + column*column_width, height_array[Index1]/2+separation*row, 0>
                finish {
                    ambient 0.1
                    diffuse 0.6
                    }
                pigment { Yellow filter 0.5 }
                //texture { PinkAlabaster }
              }
            #declare Index1 = Index1+1;
        #end
    #declare column = column + 1;
    #end
#declare row = row + 1;
#end

#declare substrate_thickness = 2;
#declare total_width = num_columns * column_width;
#declare total_height = num_rows * separation;

box{
    <-separation/2, -separation/2, thick>, <total_width-separation/2, total_height-separation/2, thick+substrate_thickness>
    //pigment { Gray filter 0.9 }
    texture { Glass }
    finish {
                    ambient 0.1
                    diffuse 2
                    }
    }
    
//cylinder{<total_width/2, total_height/2, thick  + substrate_thickness>, <total_width/2, total_height/2, thick+ substrate_thickness + 5000>, 0.9*total_height/2
//    pigment{ Green filter 0.75} 

//}

#declare height = 0.03 * total_height;
#declare outer_rad = 1.35*total_width;
#declare inner_rad = 0.95 * outer_rad;
#declare ang = 30;
#declare arrow_length = 8;
#declare arrow_width = 9;


#merge{    
    #merge{    
        #difference {
            cylinder { <0,0,0>, <0, height, 0>, outer_rad } 
            cylinder { <0,-0.01,0>, <0, 1.01*height, 0>, inner_rad }
            box { <0, -0.01, 0>, <1.001*outer_rad, 1.01*height, 1.001*outer_rad> }
            box { <0, -0.01, 0>, <-1.001*outer_rad, 1.01* height, 1.001*outer_rad> }
            box { <0, -0.01, 0>, <-1.001*outer_rad, 1.01* height, -1.001*outer_rad> }
            prism { linear_spline
             -0.01, 1.01*height, 3,
                <0, 0>, <outer_rad/tand(ang), -1.001*outer_rad>, <0, -1.001*outer_rad> }
           
            
             }
        prism { linear_spline
             -0.01, 1.01*height, 3,
                <(inner_rad+outer_rad)/2-arrow_width/2, 0>, <(inner_rad+outer_rad)/2+arrow_width/2, 0>, <(inner_rad+outer_rad)/2, -arrow_length>
                rotate ang*y
                }
        
        #declare time = 0;
        #declare aspect_ratio1 = 0.35;
        #declare alpha1 = 85;                     
        
        
    
        #merge {
            #while(time < 1)
                #declare xcoord = ellipse_rad * cos(2*pi*time);
                #declare ycoord = ellipse_rad * aspect_ratio1 * sin(2*pi*time);
                sphere {
                    <xcoord, ycoord, 0>, sphere_rad
                    
                    }
                #declare time = time + 0.001;    
            
            
            #end
            
            cone {
                <0, ellipse_rad*aspect_ratio1, 0>, arrow_rad, <arrow_height, ellipse_rad*aspect_ratio1, 0>, 0
                } 
        
            rotate alpha1*z
            rotate 90*x
            translate <(inner_rad+outer_rad)/2, height/2, 0>
            rotate ellipse_displacement*ang*y
            texture{
                pigment { Blue }
                finish { diffuse 0.75 phong 0.3}
            }
            no_shadow
        
        
        }
        
        translate <0.9*total_width/2 - outer_rad, total_height/2, move_forward>
        pigment { Green filter laser_filter }
        no_shadow
           
}
       
#merge {    
    #merge{    
        #difference {
            cylinder { <0,0,0>, <0, height, 0>, outer_rad } 
            cylinder { <0,-0.01,0>, <0, 1.01*height, 0>, inner_rad }
            box { <0, -0.01, 0>, <1.001*outer_rad, 1.01*height, 1.001*outer_rad> }
            box { <0, -0.01, 0>, <-1.001*outer_rad, 1.01* height, 1.001*outer_rad> }
            box { <0, -0.01, 0>, <-1.001*outer_rad, 1.01* height, -1.001*outer_rad> }
            prism { linear_spline
             -0.01, 1.01*height, 3,
                <0, 0>, <outer_rad/tand(ang), -1.001*outer_rad>, <0, -1.001*outer_rad> }
           
            
             }
        prism { linear_spline
             -0.01, 1.01*height, 3,
                <(inner_rad+outer_rad)/2-arrow_width/2, 0>, <(inner_rad+outer_rad)/2+arrow_width/2, 0>, <(inner_rad+outer_rad)/2, -arrow_length>
                rotate ang*y
                }
        translate <-(outer_rad+inner_rad)/2, -height/2, 0>
        rotate 180*z
        translate <(outer_rad+inner_rad)/2, height/2, 0>
        no_shadow 
           }
    #declare time = 0;
    #declare aspect_ratio2 = 0.95;
    #declare alpha2 = 40;                        
            
   #merge {
        #while(time < 1)
            #declare xcoord = ellipse_rad * cos(2*pi*time);
            #declare ycoord = ellipse_rad * aspect_ratio2 * sin(2*pi*time);
            sphere {
                <xcoord, ycoord, 0>, sphere_rad
                
                }
            #declare time = time + 0.001;    
        
        
        #end
        
        cone {
            <0, ellipse_rad*aspect_ratio2, 0>, arrow_rad, <arrow_height, ellipse_rad*aspect_ratio2, 0>, 0
            } 
    
        rotate alpha2*z
        rotate 90*x
        translate <(inner_rad+outer_rad)/2, height/2, 0>
        rotate ellipse_displacement*ang*y
        translate <-(outer_rad+inner_rad)/2, -height/2, 0>
        rotate 180*z
        translate <(outer_rad+inner_rad)/2, height/2, 0>
        texture{
                pigment { Red }
                finish { diffuse 0.75 phong 0.1}
            }
    
    
    }
    
        translate <0.9*total_width/2 - outer_rad, total_height/2, move_forward>
        pigment { Green filter laser_filter }
        no_shadow  
    
}       


    

#declare height = 0.03 * total_height;
#declare outer_rad = 0.7*total_width;
#declare inner_rad = 0.9 * outer_rad;
#declare ang = 55;
#declare arrow_length = 8;
#declare arrow_width = 9;
#merge{
    #merge{    
        #difference {
            cylinder { <0,0,0>, <0, height, 0>, outer_rad } 
            cylinder { <0,-0.01,0>, <0, 1.01*height, 0>, inner_rad }
            box { <-0.001, -0.01, 0>, <1.001*outer_rad, 1.01*height, 1.001*outer_rad> }
            box { <0.001, -0.01, 0>, <-1.001*outer_rad, 1.01* height, 1.001*outer_rad> }
            box { <0.001, -0.01, 0>, <-1.001*outer_rad, 1.01* height, -1.001*outer_rad> }
            prism { linear_spline
             -0.01, 1.01*height, 3,
                <0, 0>, <outer_rad/tand(ang), -1.001*outer_rad>, <0, -1.001*outer_rad> }
           
            
         }
        prism { linear_spline
             -0.01, 1.01*height, 3,
                <(inner_rad+outer_rad)/2-arrow_width/2, 0>, <(inner_rad+outer_rad)/2+arrow_width/2, 0>, <(inner_rad+outer_rad)/2, -arrow_length>
                rotate ang*y
         }
     }
        
        #declare time = 0;
        #declare aspect_ratio2 = 1;
        #declare alpha2 = 40*pi/180;                        
        
    
         
           
            
            
    #merge {
        #while(time < 1)
            #declare xcoord = ellipse_rad * cos(2*pi*time);
            #declare ycoord = ellipse_rad * aspect_ratio2 * sin(2*pi*time);
            sphere {
                <xcoord, ycoord, 0>, sphere_rad
                
                }
            #declare time = time + 0.001;    
        
        
        #end
        
        cone {
            <0, ellipse_rad*aspect_ratio2, 0>, arrow_rad, <arrow_height, ellipse_rad*aspect_ratio2, 0>, 0
            } 
     
        rotate alpha2*z
        rotate 90*x
        translate <(inner_rad+outer_rad)/2, height/2, 0>
        rotate ellipse_displacement*ang*y
        texture {
            pigment { Magenta }
            finish{ diffuse 0.75 phong 0.1 }
        }
    
    }
        
        translate <0.9*total_width/2 - outer_rad, total_height/2, move_forward>
        pigment { Green filter laser_filter }
        no_shadow   
}
       
#merge{    
    #difference {
         cylinder { <0,0,0>, <0, height, 0>, outer_rad } 
         cylinder { <0,-0.01,0>, <0, 1.01*height, 0>, inner_rad }
         box { <-0.001, -0.01, 0>, <1.001*outer_rad, 1.01*height, 1.001*outer_rad> }
         box { <0.001, -0.01, 0>, <-1.001*outer_rad, 1.01* height, 1.001*outer_rad> }
         box { <0.001, -0.01, 0>, <-1.001*outer_rad, 1.01* height, -1.001*outer_rad> }
        prism { linear_spline
         -0.01, 1.01*height, 3,
            <0, 0>, <outer_rad/tand(ang), -1.001*outer_rad>, <0, -1.001*outer_rad> }
       
        
     }
    prism { linear_spline
         -0.01, 1.01*height, 3,
            <(inner_rad+outer_rad)/2-arrow_width/2, 0>, <(inner_rad+outer_rad)/2+arrow_width/2, 0>, <(inner_rad+outer_rad)/2, -arrow_length>
            rotate ang*y
            }
            
    #declare time = 0;
        #declare aspect_ratio2 = 0.1;
        #declare alpha2 = 40*pi/180;                        
        
         
           
            
            
    #merge {
        #while(time < 1)
            #declare xcoord = ellipse_rad * cos(2*pi*time);
            #declare ycoord = ellipse_rad * aspect_ratio2 * sin(2*pi*time);
            sphere {
                <xcoord, ycoord, 0>, sphere_rad
                
                }
            #declare time = time + 0.001;    
        
        
        #end
        
        cone {
            <0, ellipse_rad*aspect_ratio2, 0>, arrow_rad, <arrow_height, ellipse_rad*aspect_ratio2, 0>, 0
            } 
     
        rotate alpha2*z
        rotate 90*x
        translate <(inner_rad+outer_rad)/2, height/2, 0>
        rotate ellipse_displacement*ang*y
        texture{
            pigment{ Cyan }
            finish{ diffuse 0.75 phong 0.1 }
        }
    
    
    }
    translate <-(outer_rad+inner_rad)/2, -height/2, 0>
    rotate 180*z
    translate <(outer_rad+inner_rad)/2, height/2, 0>
    translate <0.9*total_width/2 - outer_rad, total_height/2, move_forward>
    pigment { Green filter laser_filter }
    no_shadow   
       }              
 }
  light_source { <3*total_width/2, 1.25*total_height, -30> color White }

  
  camera {
    location <-20, total_height*1.35, -125>
    look_at  <total_width/2, total_height/2,  0.0>
  //  focal_point <-6, 1, 30>    // blue cylinder in focus
  //  focal_point < 0, 1,  0>    // green box in focus
    //focal_point < 1, 1, -6>    // pink sphere in focus
    aperture 0.4     // a nice compromise
  //  aperture 0.05    // almost everything is in focus
  //  aperture 1.5     // much blurring
  //  blur_samples 4       // fewer samples, faster to render
    ///blur_samples 20      // more samples, higher quality image
  }   
     
 