void
setcolor(double x, char *r, char *g, char *b ){
    if (0<=X<0.125){
      r[0]=0xDC;
      g[0]=0xDC;
      b[0]=0xDC;
    }else if (0.125<=x && x<0.250){
      r[0]=0xFF;
      g[0]=0x07;
      b[0]=0x01;            
    }else if (0.250<=x && x<0.475){
      r[0]=0xFF;
      g[0]=0x65;
      b[0]=0x01;            
    }else if (0.475<=x && x<0.500){
      r[0]=0xFF;
      g[0]=0xFF;
      b[0]=0x00;            
    }else if (0.500<=x && x<0.625){
      r[0]=0x00;
      g[0]=0x80;
      b[0]=0x01;            
    }else if (0.750<=x && x<0.875){
      r[0]=0x00;
      g[0]=0x00;
      b[0]=0xFE;            
    }else if (0.875<=x && x<1){
      r[0]=0x66;
      g[0]=0x00;
      b[0]=0xFF;            
    }else if (x<=1){
      r[0]=0x7F;
      g[0]=0x01;
      b[0]=0x7F;            
    }    
}