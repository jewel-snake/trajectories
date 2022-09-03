use std::vec::Vec;
use std::fs;
use std::ops;
use std::io::Write;

//static G: f64 = 4.30091252525e-3 * 1.022018688e-6 * 1.022018688e-6;

#[derive(Copy,Clone)]
struct Vector6d {
    x: f64,
    y: f64,
    z: f64,
    m: f64,
    n: f64,
    k: f64
}

impl Vector6d {
    fn position_mp(self, scal: f64) -> Vector6d {
        Vector6d {
            x: self.x*scal,
            y: self.y*scal,
            z: self.z*scal,
            m: self.m,
            n: self.n,
            k: self.k
        }
    }
    fn velocity_mp(self, scal: f64) -> Vector6d {
        Vector6d {
            x: self.x,
            y: self.y,
            z: self.z,
            m: self.m*scal,
            n: self.n*scal,
            k: self.k*scal
        }
    }
    fn position_div(self, scal: f64) -> Vector6d {
        if scal == 0_f64 {
            panic!("Attention! Division by 0 encountered!");
        }
        Vector6d {
            x: self.x/scal,
            y: self.y/scal,
            z: self.z/scal,
            m: self.m,
            n: self.n,
            k: self.k
        }
    }
    fn velocity_div(self, scal: f64) -> Vector6d {
        if scal == 0_f64 {
            panic!("Attention! Division by 0 encountered!");
        }
        Vector6d {
            x: self.x,
            y: self.y,
            z: self.z,
            m: self.m/scal,
            n: self.n/scal,
            k: self.k/scal
        }
    }
    fn position_sum(self,additive : Vector6d) -> Vector6d {
        Vector6d {
            x: self.x+additive.x,
            y: self.y+additive.y,
            z: self.z+additive.z,
            m: self.m,
            n: self.n,
            k: self.k
        }
    }
    fn position_sub(self, substructive: Vector6d) -> Vector6d {
        Vector6d {
            x: self.x-substructive.x,
            y: self.y-substructive.y,
            z: self.z-substructive.z,
            m: self.m,
            n: self.n,
            k: self.k
        }
    }
    fn velocity_sum(self,additive : Vector6d) -> Vector6d {
        Vector6d {
            x: self.x,
            y: self.y,
            z: self.z,
            m: self.m+additive.m,
            n: self.n+additive.n,
            k: self.k+additive.k
        }
    }
    fn velocity_sub(self, substructive: Vector6d) -> Vector6d {
        Vector6d {
            x: self.x,
            y: self.y,
            z: self.z,
            m: self.m-substructive.m,
            n: self.n-substructive.n,
            k: self.k-substructive.k
        }
    }
    fn position_mag(self) -> f64{
        (self.x*self.x+self.y*self.y+self.z*self.z).sqrt()
    }
    fn velocity_mag(self) -> f64{
        (self.m*self.m+self.n*self.n+self.k*self.k).sqrt()
    }
}
impl ops::Add<Vector6d> for Vector6d {
    type Output = Vector6d;
    fn add(self,additive : Vector6d) -> Vector6d {
        Vector6d {
            x: self.x+additive.x,
            y: self.y+additive.y,
            z: self.z+additive.z,
            m: self.m+additive.m,
            n: self.n+additive.n,
            k: self.k+additive.k
        }
    }
}

impl ops::Mul<f64> for Vector6d {
    type Output = Vector6d;
    fn mul(self, scal: f64) -> Vector6d {
        Vector6d {
            x: self.x*scal,
            y: self.y*scal,
            z: self.z*scal,
            m: self.m*scal,
            n: self.n*scal,
            k: self.k*scal
        }
    }
}

impl ops::Div<f64> for Vector6d {
    type Output = Vector6d;
    fn div(self, scal: f64) -> Vector6d {
        if scal == 0_f64 {
            panic!("Attention! Division by 0 encountered!");
        }
        Vector6d {
            x: self.x/scal,
            y: self.y/scal,
            z: self.z/scal,
            m: self.m/scal,
            n: self.n/scal,
            k: self.k/scal
        }
    }
}

impl ops::Sub<Vector6d> for Vector6d {
    type Output = Vector6d;
    fn sub(self, substructive: Vector6d) -> Vector6d {
        Vector6d {
            x: self.x-substructive.x,
            y: self.y-substructive.y,
            z: self.z-substructive.z,
            m: self.m-substructive.m,
            n: self.n-substructive.n,
            k: self.k-substructive.k
        }
    }
}

/*fn law(r: f64) -> f64 {
    -G*1e12/(17500_f64).powi(2) -1.364289272
}*/

fn law(instance: Vector6d) -> f64 {
    let r = (instance.x.powi(2)+instance.y.powi(2)).sqrt();
    fi_b(instance.position_mag()) + fi_d(r,instance.z) + fi_h(instance.position_mag())
}

fn fi_b(r: f64) -> f64{
    -142.0/(r.powi(2)+249.7_f64.powi(2)).sqrt()
}

fn fi_d(r: f64,z:f64) -> f64{
    -2732.0/(r.powi(2)+(5160.0+(z.powi(2)+310.5_f64.powi(2)).sqrt()).powi(2)).sqrt()
}

fn fi_h(r: f64) -> f64{
    -24572.0/64300_f64*((64300_f64+(64300_f64.powi(2)+r.powi(2)).sqrt())/r).ln()
}

fn dm_dt(instance: Vector6d) -> Vector6d {
    let r = instance.position_mag();
    let accel = law(instance) * 1.0749776097859195e-11;// * 1.022018688e-6_f64.powi(2)*2.325e7_f64/100.0;
    Vector6d {
        x: instance.m,
        y: instance.n,
        z: instance.k,
        m: accel * (instance.x/r),
        n: accel * (instance.y/r),
        k: accel * (instance.z/r)
    }
}

fn integrate(start: Vector6d, data: &mut Vec<Vector6d>, whole_time: f64) /*-> Result<Vec<Vector6d>,str>*/{
    //let mut t = 0;
    //let mut data = Vec::with_capasity(i);
    let i = data.capacity();
    let incr = whole_time / i as f64;
    data.push(start);
    for j in 1..i {
        let m = data[j-1];
        let k1 = dm_dt(m);
        let k2 = dm_dt(m+k1*(incr/2_f64));
        let k3 = dm_dt(m+k2*(incr/2_f64));
        let k4 = dm_dt(m+k3*incr);
        data.push(m+(k1+k2*2.0+k3*2.0+k4)*(incr/6.0));
    }
    //Ok(data)
}

fn main() {
    let start_data = fs::read_to_string("./initdata").unwrap().trim().split(',').map(|x| x.parse::<f64>().unwrap()).collect::<Vec<f64>>();
    let start = Vector6d {
        x: start_data[0] * 1000_f64,
        y: start_data[1] * 1000_f64,
        z: start_data[2] * 1000_f64,
        m: start_data[3] * 1.022018688e-6, //1.0219053791315e-6,
        n: start_data[4] * 1.022018688e-6, //1.0219053791315e-6,
        k: start_data[5] * 1.022018688e-6  //1.0219053791315e-6
    };
    let mut trace = Vec::with_capacity(10000);
    integrate(start,&mut trace,2e9);
    let name = &std::env::args().collect::<Vec<String>>()[1];
    let mut file = fs::File::create(name).unwrap();
    for i in trace {
        file.write(format!("{:.9},{:.9},{:.9}\n", i.x/1000.0,i.y/1000.0,i.z/1000.0).as_bytes());
    }
}
