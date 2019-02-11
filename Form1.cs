using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Net;
namespace _3dpixels
{
    
    public partial class Form1 : Form
    {
        public Bitmap bm;
        public Graphics graphics;
        public Camera Camera1;
        ParticleSystem PSystem;
        public Point CameraAngles;
        public Plane PlaneA;
        public int PreviousFrame = 0;
        public List<Mesh> Cubes;
        public LightSource Sun;

        public Mesh Sphere;
       
        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            Sphere = new Mesh();
            Sphere.ToUVSphere(2, 30);
            Camera1 = new Camera();
            Camera1.FOV = new Size(90, 90);
            Camera1.NormalAxis = Axes.Z;
            Camera1.VerticalAxis = Axes.Y;
            Camera1.RightAxis = Axes.X;
            Cubes = new List<Mesh> { new Mesh(), new Mesh(), new Mesh() };
            Sun = new LightSource();
            Sun.Position = new _Point(0, 10, 1);
            bm = new Bitmap(800, 800);
            Globals.BMSize = bm.Size;
            graphics = Graphics.FromImage(bm);
            for (int x=0; x<bm.Width; x++)
            {
                for(int y=0; y<bm.Height; y++)
                {
                    bm.SetPixel(x, y, Color.White);
                }
            }
         

            CameraAngles = new Point(90, 90);

            /*
            PSystem = new ParticleSystem();

            PSystem.GenerateParticles(10, new _Point(0, 0, 4), 1, new _Point(0, 0, 0));
            PSystem.TimeInterval = .1;
            PSystem.ParticleRadius = .01;
            PSystem.AttractionRadius = .5;
            PSystem.ParticleColor = Color.Black;
            DrawParticles(PSystem);
            */
            Cubes[0].ToCube(new _Point(-1, -1, 2), .8);
            Cubes[1].ToCube(new _Point(0, 1.5, 2), 1.5);
            Cubes[2].ToCube(new _Point(2, -1, 3), 1);

              foreach (Mesh cube in Cubes)
             {
                  DrawMesh(cube);
               }


            //DrawMesh(Sphere);
            /*
            for(int i=0; i<Sphere.Verts.Count; i++)
            {
                for (int j = 0; j < Sphere.Verts[i].Count; j++)
                {
                    _Point point = Sphere.Verts[i][j];
                    Sphere.Verts[i][j] = VectMath.Add(new _Point(0, -2, 4), point);
                    point = Sphere.Verts[i][j];
                    if (convert2D(point).X > 0 && convert2D(point).X < bm.Width && convert2D(point).Y > 0 && convert2D(point).X < bm.Height)
                    {

                        bm.SetPixel(convert2D(point).X, convert2D(point).Y, Color.Black);
                    }
                }
            }
            */
            pictureBox1.Image = bm;
        }

        public void DrawPlane(Plane plane)
        {
            if (VectMath.AngleBetween(VectMath.Subtract(plane.Position, new _Point(0, 0, 0)), plane.Normal) > Math.PI/2)
            {
                PointF[] pf = new PointF[plane.Points.Count];

                for (int i = 0; i < pf.Length; i++)
                {
                    Point p = convert2D(plane.Points[i]);
                    pf[i] = new PointF(p.X, p.Y);

                }
                graphics.FillPolygon(plane.PerceivedColor(Sun, new _Point(0, 0, 0)), pf);
                pictureBox1.Image = bm;
               // Globals.PlanePointsToClear.Add(pf);
            }
           // MessageBox.Show(VectMath.AngleBetween(VectMath.Subtract(plane.Position, new _Point(0, 0, 0)), plane.Normal).ToString());
        }
        public void DrawMesh(Mesh mesh)
        {
            foreach(Plane plane in mesh.Planes)
            {
                
                DrawPlane(plane);
            }
        }
        public void DrawParticles(ParticleSystem ps)
        {
            
            foreach(Particle particle in ps.Particles)
            {
                Point p = new Point();
               
                p = convert2D(particle.Position);               
                bm.SetPixel(p.X, p.Y, ps.ParticleColor);
                
            }
            pictureBox1.Image = bm;
        }

        public void ClearBM()
        {
            graphics.Clear(Color.Black);

        }

        private void button1_Click(object sender, EventArgs e)
        {
            Mesh cube = new Mesh();
            cube.ToCube(new _Point(8*Globals.rnd.NextDouble()-4, 8 * Globals.rnd.NextDouble() - 4, 14),3*Globals.rnd.NextDouble());
            Cubes.Add(cube);

            
        }

        public void button2_Click(object sender, EventArgs e)
        {
            graphics.Clear(Color.White);



            foreach (Mesh cube in Cubes)
             {
                cube.Rotate(new _Point(Convert.ToInt32(textBox4.Text), Convert.ToInt32(textBox5.Text), Convert.ToInt32(textBox6.Text)), 10);

               DrawMesh(cube);
             }
            Sphere.Rotate(new _Point(Convert.ToInt32(textBox4.Text), Convert.ToInt32(textBox5.Text), Convert.ToInt32(textBox6.Text)), 10);
            DrawMesh(Sphere);

            pictureBox1.Image = bm;

            //PlaneA.Rotate(Axes.X, 30);
            
            //DrawPlane(PlaneA);
            
        }

        private void button3_Click(object sender, EventArgs e)
        {
            _Point p = new _Point(Convert.ToDouble(textBox1.Text), Convert.ToDouble(textBox2.Text), Convert.ToDouble(textBox3.Text));
            _Point v = new _Point(Convert.ToDouble(textBox4.Text), Convert.ToDouble(textBox5.Text), Convert.ToDouble(textBox6.Text));
            label4.Text = VectMath.RotateByAngle(p, v, Convert.ToDouble(textBox7.Text)).x.ToString() + ", " + VectMath.RotateByAngle(p, v, Convert.ToDouble(textBox7.Text)).y.ToString() + ", " + VectMath.RotateByAngle(p, v, Convert.ToDouble(textBox7.Text)).z.ToString();


        }
        int frame = 0;
        private void timer1_Tick(object sender, EventArgs e)
        {
            frame += 1;
            graphics.Clear(Color.White);
            //Sphere.Rotate(new _Point(0, 1, 0), 5);
            //Sphere.SetPosition(VectMath.Add(Sphere.Position, VectMath.Scale(new _Point(0, 1, 0), .02*Math.Sin(frame*.1))));
            //DrawMesh(Sphere);
            Cubes[0].Rotate(Axes.Y, 6);
            Cubes[1].Rotate(Axes.X, 6);
            Cubes[2].Rotate(Axes.Z, 6);
            
            
            foreach (Mesh cube in Cubes)
            {
                if (cube.Position.z > 1)
                {
                    cube.SetPosition(VectMath.Add(cube.Position, new _Point(0, 0, -.1*Math.Sqrt(Math.Pow(cube.Position.x, 2) + Math.Pow(cube.Position.y,2)))));
                }
                else { cube.SetPosition(new _Point(8 * Globals.rnd.NextDouble()-4, 8 * Globals.rnd.NextDouble()-4, 15)); }
                DrawMesh(cube);
                
            }
            
        }
        public Point convert2D(_Point p)
        {
            if(Camera1.NormalAxis != Axes.Z)
            {
                //IO.Point(p);
                p = VectMath.RotateByAngle(p, VectMath.Cross(Axes.Z, Camera1.NormalAxis), 180*(VectMath.AngleBetween(new _Point(0, 0, 1), Camera1.NormalAxis))/Math.PI);
                //IO.Point(p);
            }
             //if( Camera1.VerticalAxis != Axes.Y)
           // {
           //     p = VectMath.RotateByAngle(p, Camera1.NormalAxis, 180*(VectMath.AngleBetween(new _Point(0, 1, 0), Camera1.VerticalAxis))/Math.PI);
           // }
             if(p.z<=0)
            {
                return (new Point(0, 0));
            }
            int Xcoord = Convert.ToInt32((p.x * Globals.BMSize.Width / 2) / (p.z * Math.Tan(45))) + Globals.BMSize.Width / 2;
            int Ycoord = Convert.ToInt32((-p.y * Globals.BMSize.Height / 2) / (p.z * Math.Tan(45))) + Globals.BMSize.Height / 2;
            p.PF = new PointF(Xcoord, Ycoord);
            
            return (new Point(Xcoord, Ycoord));
        }

        private void pictureBox2_Click(object sender, EventArgs e)
        {
            RotateTimer.Enabled = true;

            //Sphere.Rotate(new _Point(0, 1, 0), 5);
            //Sphere.SetPosition(VectMath.Add(Sphere.Position, VectMath.Scale(new _Point(0, 1, 0), .02*Math.Sin(frame*.1))));
            //DrawMesh(Sphere);
            //Cubes[0].Rotate(Axes.Y, 6);
            //Cubes[1].Rotate(Axes.X, 6);
           // Cubes[2].Rotate(Axes.Z, 6);


           // foreach (Mesh cube in Cubes)
           // {
           //     if (cube.Position.z > 1)
           //     {
           //         cube.SetPosition(VectMath.Add(cube.Position, new _Point(0, 0, -.1 * Math.Sqrt(Math.Pow(cube.Position.x, 2) + Math.Pow(cube.Position.y, 2)))));
           //     }
          //      else { cube.SetPosition(new _Point(8 * Globals.rnd.NextDouble(), 8 * Globals.rnd.NextDouble(), 15)); }
           //     DrawMesh(cube);

         //   }
        }

        private void timer2_Tick(object sender, EventArgs e)
        {
            Point Vmouse = new Point(MousePosition.X - Globals.LastMPosition.X, MousePosition.Y - Globals.LastMPosition.Y);
            if(Math.Sqrt(Vmouse.X^2+Vmouse.Y^2)>50)
            {

                //Camera1.Rotate()
            }
            Globals.LastMPosition = MousePosition;
        }

        private void button4_Click(object sender, EventArgs e)
        {
            Globals.InterfaceMode = "Edit";
            
        }

        private void pictureBox1_Click(object sender, EventArgs e)
        {
            if(Globals.InterfaceMode == "Edit")
            {
                Point MouseLoc = System.Windows.Forms.Cursor.Position;
                int numpts = 0;
                foreach(PointF p in Globals.PointFs)
                {

                }
            }
        }

        private void RotateTimer_Tick(object sender, EventArgs e)
        {
            Globals.Frame += 1;
            Camera1.Rotate(Axes.X, 5);
            
            graphics.Clear(Color.White);
            foreach (Mesh cube in Cubes)
            {
                DrawMesh(cube);
            }
        }

        private void button5_Click(object sender, EventArgs e)
        {
            FPSTimer.Enabled = true;
        }

        private void FPSTimer_Tick(object sender, EventArgs e)
        {
            int FPS = Globals.Frame - PreviousFrame;
            PreviousFrame = Globals.Frame;
            label6.Text = FPS.ToString();
        }
    }
}
public class _Point
{
    public double x;
    public double y;
    public double z;
    public PointF PF;
    public _Point(double x, double y, double z)
    {
        this.x = x;
        this.y = y;
        this.z = z;
        
    }
    
}
public class Plane
{
    public List<_Point> Points;
    public _Point Normal;
    public _Point LightVector;
    public Color Color;
    public _Point Position;
    
    public Plane(List<_Point> Points)
    {
        this.Points = Points;
        Position = VectMath.CenterofPoints(Points);
        Normal = CalculateNormal();
    }
    public void SetPosition(_Point p)
    {
        //IO.Point(Position);
        _Point DiffVector = VectMath.Subtract(p, Position);
        for(int i=0; i < Points.Count; i++)
        {
            Points[i] = VectMath.Add(Points[i], DiffVector);
        }
        Position = p;
        Normal = CalculateNormal();
    }
    public void Rotate(_Point axis, double deg)
    {
        for(int i=0; i<Points.Count; i++)
        {
            Points[i] = VectMath.Add(VectMath.RotateByAngle(VectMath.Subtract(Points[i], Position), axis, deg), Position);
        }
        Normal = CalculateNormal();
    }
    public Brush PerceivedColor(LightSource light, _Point CamPosition)
    {
        
        _Point LightVector = VectMath.Subtract(this.Position, light.Position);
        _Point PercVector = this.Position;
        _Point Rvector = VectMath.RotateByAngle(VectMath.Negate(LightVector), Normal, 180);
        double whitefactor = Math.Abs(.5*Math.Sin(VectMath.AngleBetween(Normal, LightVector))+.5);
        double colorfactor = Math.Abs( Math.Cos(VectMath.AngleBetween(Normal, LightVector)));
        double Total = Math.Pow(Math.Sin(VectMath.AngleBetween(PercVector, Rvector)), 2);
        SolidBrush brush = new SolidBrush( Color.FromArgb(Convert.ToInt32(255 * Total * colorfactor * whitefactor), Convert.ToInt32(255 * Total * whitefactor), Convert.ToInt32(255 * Total * whitefactor)));
        return (brush);
    }
    public _Point CalculateNormal()
    {
        return(VectMath.Cross(VectMath.Subtract(Points[1], Points[0]), VectMath.Subtract(Points[2], Points[0])));
    }
}

public class ParticleSystem
{
    public double ParticleRadius;
    public double AttractionRadius;
    public double TimeInterval;
    
    public Color ParticleColor;
    public List<Particle> Particles;
    public void GenerateParticles(int num, _Point origin, double radius, _Point velocity)
    {
        Particles = new List<Particle>();
        for(int i=0; i<num; i++)
        {
            Particles.Add(new Particle(new _Point(origin.x+ Globals.rnd.NextDouble()*radius, origin.y+Globals.rnd.NextDouble() * radius, origin.z + Globals.rnd.NextDouble() * radius)));
            //MessageBox.Show(Particles[i].Position.x.ToString());
            Particles[i].velocity = velocity;
            Particles[i].force = new _Point(0, 0, 0);
            Particles[i].Mass = 1;
        }
    }
    public void Iterrate()
    {
        foreach(Particle particle in Particles)
        {
            particle.force = new _Point(0, 0, 0);
            foreach(Particle p in Particles)
            {
                particle.force = VectMath.Add(ForceOn(particle, p), particle.force);
            }
            particle.velocity = VectMath.Add(particle.velocity, (VectMath.Scale(VectMath.Scale(particle.force, 1 / particle.Mass), TimeInterval)));
            //MessageBox.Show(VectMath.Magnitude(particle.force).ToString());
            
        }
        foreach(Particle particle in Particles)
        {
            particle.Position = VectMath.Add(particle.Position, VectMath.Scale(particle.velocity, TimeInterval));
        }
    }
    public _Point ForceOn(Particle particle, Particle p)
    {
        _Point vector = VectMath.Subtract(p.Position, particle.Position);      //Vector pointing from particle to p
        double distance = VectMath.Magnitude(vector);
        if(distance > ParticleRadius & distance < AttractionRadius)
        {
            
            return (VectMath.Scale(vector, (1 / distance)));
                       
        }
        else
        {
            if(distance < ParticleRadius)
            {
                if (distance == 0)
                {
                    return (new _Point(0, 0, 0));
                }
                else { return (VectMath.Scale(vector, (-1 / distance))); }
            }
            else { return (new _Point(0,0,0)); }
        }
    

    }


}

public class Particle
{
    public _Point force;
    public _Point velocity;
    public _Point Position;
    public double Mass;
    public Particle(_Point p)
    {
        Position = p;
    }
}

public class Globals
{
    public static Random rnd = new Random();
    public static List<_Point> Points3D = new List<_Point>();
    public static List<PointF> PointFs = new List<PointF>();
    public static Size BMSize = new Size(900, 900);
    public static Point LastMPosition = new Point(0, 0);
    public static int Frame = 0;
    public static int NumFramesLastSecond = 0;
    public static String InterfaceMode = "";

}

public class Mesh
{
    public List<List<_Point>> Verts;
    public List<Plane> Planes;
    public _Point Position;
    public void SetPosition(_Point p)
    {
        _Point MoveVector = VectMath.Subtract(p, Position);
        foreach(Plane plane in Planes)
        {
            plane.SetPosition(VectMath.Add(plane.Position, MoveVector));
        }
        Position = p;
    }
    public void ToCube(_Point position, double l)
    {
        this.Position = new _Point(0,0,l/2);
        if (Planes != null)
        {
            Planes.Clear();
        }
        else
        {
            Planes = new List<Plane>();
        }
                
        for (int i = 0; i < 6; i++)
        {
            Planes.Add(new Plane(new List<_Point> { new _Point(0,0,0), new _Point(0,l,0), new _Point(l,l,0), new _Point(l,0,0)}));
        }
        Planes[0].SetPosition(new _Point(0, 0, 0));
        //Planes[0].Rotate(Axes.Y, 180);
        Planes[1].SetPosition(new _Point(0, l / 2, l / 2));
       
        Planes[1].Rotate(Axes.X, 90);
        Planes[2].SetPosition(new _Point(0, 0, l));
        Planes[2].Rotate(Axes.Y, 180);
        Planes[3].SetPosition(new _Point(0, -l / 2, l / 2));
        Planes[3].Rotate(Axes.X, 270);
        Planes[4].SetPosition(new _Point(-l / 2, 0, l / 2));
        Planes[4].Rotate(Axes.Y, 90);
        Planes[5].SetPosition(new _Point(l / 2, 0, l / 2));
        Planes[5].Rotate(Axes.Y, 270);

        this.SetPosition(position);
    }

    public void ToUVSphere(double radius, int Numrings)
    {
        Verts = new List<List<_Point>>();
        Position = new _Point(0, 0, 0);
        Planes = new List<Plane>();
        double angle = (180 / Numrings)*(Math.PI/180);
        double VertAngle = (180 / Numrings) * (Math.PI / 180);
        //for (int i=0; i<2*Numrings; i++)
        //{
        //    Plane plane = new Plane(new List<_Point> {new _Point(0,radius,0), new _Point(radius * Math.Cos(VertAngle) * Math.Sin(angle * (i+1)), radius * Math.Cos(VertAngle), radius * Math.Cos(VertAngle) * -Math.Cos(angle * (i+1))), new _Point(radius * Math.Cos(VertAngle) * Math.Sin(angle * i), radius * Math.Cos(VertAngle), radius * Math.Cos(VertAngle) * -Math.Cos(angle * i)) });
        //    Planes.Add(plane);
       // }
        
        for (int i = 0; i <  Numrings; i++)
        {
            List<_Point> list = new List<_Point>();
            Verts.Add(list);
            for (int j=0; j<2*Numrings; j++)
            {
                Verts[i].Add(new _Point(radius * Math.Sin(VertAngle*i) * Math.Sin(angle * j), radius * Math.Cos(VertAngle*i), radius * Math.Sin(VertAngle*i) * -Math.Cos(angle * j)));
            }

        }
        for (int i = 0; i < Verts.Count- 1; i++)
        {
            for (int j = 0; j < Verts[i].Count-1; j++)
            {
                Plane plane = new Plane(new List<_Point> {Verts[i + 1][j], Verts[i][j], Verts[i][j+1], Verts[i+1][j+1]});
                Planes.Add(plane);
            }
            Plane p= new Plane(new List<_Point> { Verts[i+1][Verts[i].Count - 1], Verts[i][Verts[i].Count-1], Verts[i][0], Verts[i+1][0] });
            Planes.Add(p);
        }
        for(int i=0; i<Verts.Count-1; i++)
        {
            Plane p = new Plane(new List<_Point> { Verts[1][i], Verts[0][i], Verts[0][i+1], Verts[1][i+1] });
            Planes.Add(p);
            p= new Plane(new List<_Point> { Verts[Verts.Count - 1][i], Verts[Verts.Count - 2][i], Verts[Verts.Count-2][i + 1], Verts[Verts.Count - 1][i + 1] });
            Planes.Add(p);
        }
        Planes.Add(new Plane(new List<_Point> { Verts[1][Verts[1].Count - 1], Verts[0][Verts[0].Count - 1], Verts[0][0], Verts[1][0] }));

        SetPosition(new _Point(0, -1, 5));
    }



    public void Rotate(_Point axis, double angle)
    {
        foreach(Plane plane in Planes)
        {
            for(int i=0; i<plane.Points.Count; i++)
            {
                _Point p = plane.Points[i];
                //IO.Point(p);
                plane.Points[i] = VectMath.Add(VectMath.RotateByAngle(VectMath.Subtract(p, this.Position), axis, angle), this.Position);
                //IO.Point(plane.Points[i]);
            }
            plane.Position = VectMath.CenterofPoints(plane.Points);
            plane.Normal = plane.CalculateNormal();
        }
    }
    public void Move(_Point point)
    {
        SetPosition(VectMath.Add(Position, point));
    }





}


public class IO
{
    public static void Point(_Point p)
    {
        MessageBox.Show("(" + p.x.ToString() + ", " + p.y.ToString() + ", " + p.z.ToString() + ")");
    }
}

public class VectMath
{
    public static _Point Subtract(_Point a, _Point b)
    {
        return (new _Point(a.x - b.x, a.y - b.y, a.z - b.z));
    }
    public static _Point Add(_Point a, _Point b)
    {
        return (new _Point(a.x + b.x, a.y + b.y, a.z + b.z));
    }
    public static double Magnitude(_Point v)
    {
        return (Math.Sqrt(Math.Pow(v.x, 2) + Math.Pow(v.y, 2) + Math.Pow(v.z, 2)));
    }
    public static double SqrMagnitude(_Point v)
    {
        return ((Math.Pow(v.x, 2) + Math.Pow(v.y, 2) + Math.Pow(v.z, 2)));
    }
    public static _Point Scale(_Point p, double a)
    {
        return (new _Point(p.x * a, p.y * a, p.z * a));
    }
    public static double Dot(_Point a, _Point b)
    {
        return (a.x * b.x + a.y * b.y + a.z * b.z);
    }
    public static _Point Project(_Point a, _Point b)
    {
        return (VectMath.Scale(b, (VectMath.Dot(a, b) / (VectMath.SqrMagnitude(b)))));
    }
    public static _Point Cross(_Point a, _Point b)
    {
        return (new _Point(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x));
    }
    public static double AngleBetween(_Point a, _Point b)
    {
        return (Math.Acos((VectMath.Dot(a, b)) / (VectMath.Magnitude(a) * VectMath.Magnitude(b))));
    }
    public static _Point RotateByAngle(_Point pnt, _Point Axis, double Angle)
    {        
        if(VectMath.IsParallel(pnt, Axis)==true)
        {
            return (pnt);
        }
        Angle = (Angle * Math.PI) / 180;
        _Point Avector = VectMath.Subtract(pnt, VectMath.Project(pnt, Axis));
        double MagAvector = VectMath.Magnitude(Avector);
        _Point V1 = VectMath.Unit(VectMath.Cross(Axis, Avector));     
        _Point V2 = (VectMath.Negate(Avector));        
        _Point TransVector = VectMath.Add(VectMath.Scale(V1, MagAvector * Math.Sin(Angle)), VectMath.Scale(V2, (-1*Math.Cos(Angle)+1)));        
        return (VectMath.Add(TransVector, pnt));        
    }
    public static _Point Negate(_Point a)
    {
        return new _Point(-a.x, -a.y, -a.z);
    }
    
    public static _Point Unit(_Point a)
    {
        return (VectMath.Scale(a, 1 / VectMath.Magnitude(a)));
    }
    public static Boolean IsParallel(_Point a, _Point b)
    {
        if(a.x/b.x == a.y/b.y && a.y/b.y == a.z/b.z)
        {
            return (true);
        }
        else { return (false); }
    }
    public static _Point CenterofPoints(List<_Point> points)
    {
        double Netx=0;
        double Nety=0;
        double Netz=0;
        foreach(_Point p in points)
        {
            Netx += p.x;
            Nety += p.y;
            Netz += p.z;
        }
        return (new _Point((Netx / points.Count), (Nety / points.Count), (Netz / points.Count)));
    }


}

public class LightSource
{
    public string Type;
    public _Point Position;
    public Boolean IsActive;
}

public class Axes
{
    public static _Point X = new _Point(1, 0, 0);
    public static _Point Y = new _Point(0, 1, 0);
    public static _Point Z = new _Point(0, 0, 1);
}
public class Camera
{
    public _Point NormalAxis;
    public _Point VerticalAxis;
    public _Point RightAxis;
    public Size FOV;
    public Boolean Active;

    public void Rotate(_Point Axis, double Angle)
    {
        NormalAxis = VectMath.RotateByAngle(NormalAxis, Axis, Angle);
       VerticalAxis = VectMath.RotateByAngle(VerticalAxis, Axis, Angle);
    }
}