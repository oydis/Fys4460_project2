#include "Integrator.h"

Integrator::Integrator(double dt, double m, int numpart, int dim, Atom *atoms, double b, int Nc, double sig, double eps, double T_bath):
    dt(dt),
    m(m),
    numpart(numpart),
    dim(dim),
    atoms(atoms),
    b(b),
    Nc(Nc),
    sig(sig),
    eps(eps),
    T_bath(T_bath)
{
}

void Integrator::integrate()
{
    int i,t,j,k,l,g,h,o,p,q,i2,j2,k2;
    int n = 20000;
    vec v (dim);
    vec r (dim);
    vec r2 (dim);
    vec d (dim);
    vec F (dim);
    vec F_bidrag (dim);
    vec F_old (dim);
    double pot = 0.0;
    double E_k = 0.0;
    double E_tot;
    double T;
    vec v_half (dim);
    vec r_new (dim);
    vec v_new (dim);
    vec r_min (dim);
    double d_length;
    double L = b*(Nc);
    double V = L*L*L;
    double Fr_sum = 0.0;
    double P;
    vec r_init (dim);
    vec r_disp (dim);
    vec r_disp_min (dim);
    double r_disp_2;
    double r_sum = 0.0;
    double mean_disp;
    double D;
    int bins = 50;
    vec g_r (bins);
    int dist;
    double volume;
    double pi = 4*atan(1.0);
    double tau = 10*dt;
    double gamma;
    double tall;
    double T_thermalize = 1.05;
    //double T_thermalize = 0.851;
    double T_thermalize2 = 1.5;
//    double R = 4.0;
//    double R = 5.874;
//    double R = 7.0;
//    double R = 9.0;
//    double R = 11.0;
    double R = 15.0;
    bool can_move;
    bool can_move2;
    int moving_atoms = 0;
    int moving_atoms1 = 0;
    int moving_atoms2 = 0;
    int moving_atomsP1 = 0;
    int moving_atomsP2 = 0;
    int non_moving_atoms = 0;
    int spheres = 20;
    double sphere_volume;
    vec sphere_pos_x (spheres);
    vec sphere_pos_y (spheres);
    vec sphere_pos_z (spheres);
    vec sphere_radius (spheres);
    bool removed_atom;
    bool removed_atom2;
    int removed_atoms = 0;
    double V_cell;
    double E_k1 = 0;
    double E_k2 = 0;
    double P1,P2,T1,T2;
    double Fr_sum1 = 0;
    double Fr_sum2 = 0;
    double volume_dP;
    double U = 0;
    int n_u = 40;
    vec u (n_u);
    vec uz (n_u);
    vec number_of_elements (n_u);
    vec radius (n_u+1);
    vec volume_cyl (n_u);
    vec density (n_u);
    double v_xyz,r_xy;

    for(i=0;i<n_u+1;i++){
        radius(i) = i*R/30.0;
    }
    for(i=0;i<n_u;i++){
        volume_cyl(i) = pi*L*(pow(radius(i+1),2) - pow(radius(i),2));
    }

    for(int i=0;i<spheres;i++){
        sphere_radius(i) = 8.811 - 2.937*randu();
        sphere_pos_x(i) = L*randu();
        sphere_pos_y(i) = L*randu();
        sphere_pos_z(i) = L*randu();
    }

    double r_cut = 3.0;
    int cells_x = L/r_cut;
    int nx,ny,nz;
    std::vector<Atom*> celle;
    std::vector<Atom*> celle2;

    Cell_container Con(cells_x);


    for(i=0;i<numpart;i++){
        r = atoms[i].getPosition();
        nx = int (r(0)/r_cut);
        ny = int (r(1)/r_cut);
        nz = int (r(2)/r_cut);

        if(nx>cells_x-1){
           nx = cells_x-1;
        }
        if(ny>cells_x-1){
           ny = cells_x-1;
        }
        if(nz>cells_x-1){
           nz = cells_x-1;
        }

        Con.container(nx,ny,nz,&atoms[i]);
    }



    ofstream myfile;
    myfile.open ("../../results/fluid_cylinder_u_uz_rho_U_test.xyz");

    myfile << "Velocity, density, flow.\n";


    for(i=0;i<cells_x;i++){  // Finn kraft første gang
        for(j=0;j<cells_x;j++){
            for(k=0;k<cells_x;k++){
                celle = Con.getContainer(i,j,k);

                for(l=0;l<celle.size();l++){
                    r = celle[l]->getPosition();

                    for(g=0;g<dim;g++){
                        F(g) = 0;
                    }

                    for(g=0;g<celle.size();g++){  // egen celle
                        if(g!=l){
                            r2 = celle[g]->getPosition();
                            d = r - r2;
                            d_length = sqrt(d(0)*d(0) + d(1)*d(1) + d(2)*d(2));

                            pot = pot + 4*(pow (1/d_length,12) - pow (1/d_length,6));

                            for(h=0;h<dim;h++){
                                F_bidrag(h) = 24*(2*pow (1/d_length,12) - pow (1/d_length,6))*d(h)/(d_length*d_length);
                                F(h) = F(h) + F_bidrag(h);
                                Fr_sum = Fr_sum + 0.5*F_bidrag(h)*d(h);
                            }
                        }
                    }
                    for(g=-1;g<2;g++){   // naboceller
                        for(h=-1;h<2;h++){
                            for(o=-1;o<2;o++){
                                if(!(g==0 && h==0 && o==0)){

                                    i2 = i+g;
                                    j2 = j+h;
                                    k2 = k+o;

                                    if(i+g < 0){
                                       i2 = cells_x-1;
                                    }
                                    if(j+h < 0){
                                       j2 = cells_x-1;
                                    }
                                    if(k+o < 0){
                                       k2 = cells_x-1;
                                    }

                                    if(i+g == cells_x){
                                       i2 = 0;
                                    }
                                    if(j+h == cells_x){
                                       j2 = 0;
                                    }
                                    if(k+o == cells_x){
                                       k2 = 0;
                                    }

                                    celle2 = Con.getContainer(i2,j2,k2);

                                    for(p=0;p<celle2.size();p++){
                                        r2 = celle2[p]->getPosition();

                                        d = r-r2;

                                        for(q=0;q<dim;q++){   // minimum image convention
                                            if(fabs(d(q)-L) < fabs(d(q))){
                                               if(fabs(d(q)-L) < fabs(d(q)+L)){
                                                  r_min(q) = d(q)-L;
                                               }
                                            }
                                            if(fabs(d(q)) < fabs(d(q)+L)){
                                               if(fabs(d(q)) < fabs(d(q)-L)){
                                                  r_min(q) = d(q);
                                               }
                                            }
                                            if(fabs(d(q)+L) < fabs(d(q)-L)){
                                               if(fabs(d(q)+L) < fabs(d(q))){
                                                  r_min(q) = d(q)+L;
                                               }
                                            }
                                        }
                                        d_length = sqrt(r_min(0)*r_min(0) + r_min(1)*r_min(1) + r_min(2)*r_min(2));

                                        pot = pot + 4*(pow (1/d_length,12) - pow (1/d_length,6));

                                        for(q=0;q<dim;q++){
                                            F_bidrag(q) = 24*(2*pow (1/d_length,12) - pow (1/d_length,6))*r_min(q)/(d_length*d_length);
                                            F(q) = F(q) + F_bidrag(q);
                                            Fr_sum = Fr_sum + 0.5*F_bidrag(q)*r_min(q);

                                        }
                                    }
                                }
                            }
                        }
                    }
                    celle[l]->setForce(F);

                }
            }
        }
    }

    for(t=0;t<n;t++){  // integrerer


        cout << t << endl;

//        if(t>202){
//            myfile << moving_atoms << endl;
//            myfile << "Argon atoms in fcc lattice.\n";
//        }

        if(t==203){
            myfile << moving_atoms << " " << R << " " << L << "\n";
        }


        for(i=0; i<numpart;i++){  // Finner nye posisjoner

            r = atoms[i].getPosition();
            v = atoms[i].getVelocity();
            F = atoms[i].getForce();
            can_move = atoms[i].getMove();
            removed_atom = atoms[i].testRemoved();

            if(removed_atom==false && can_move==true){
                E_k = E_k + 1.0/2.0*(v(0)*v(0) + v(1)*v(1) + v(2)*v(2));
            }

//            if(can_move==true && removed_atom==false){   // Finding mean square displacement
//                r_init = atoms[i].getPosInit();

//                for(j=0;j<dim;j++){
//                    r_disp(j) = r(j) - r_init(j);
//                }

//                for(q=0;q<dim;q++){   // minimum image convention
//                    if(fabs(r_disp(q)-L) < fabs(r_disp(q))){
//                        if(fabs(r_disp(q)-L) < fabs(r_disp(q)+L)){
//                            r_disp_min(q) = r_disp(q)-L;
//                        }
//                    }
//                    if(fabs(r_disp(q)) < fabs(r_disp(q)+L)){
//                        if(fabs(r_disp(q)) < fabs(r_disp(q)-L)){
//                            r_disp_min(q) = r_disp(q);
//                        }
//                    }
//                    if(fabs(r_disp(q)+L) < fabs(r_disp(q)-L)){
//                        if(fabs(r_disp(q)+L) < fabs(r_disp(q))){
//                            r_disp_min(q) = r_disp(q)+L;
//                        }
//                    }
//                }

//                r_disp_2 = r_disp_min(0)*r_disp_min(0) + r_disp_min(1)*r_disp_min(1) + r_disp_min(2)*r_disp_min(2);

//                r_sum = r_sum + r_disp_2;
//            }


//            if(t>202 && can_move==true && removed_atom==false){
//                myfile << "Ar";

//                for(j=0;j<dim;j++){
//                    myfile << " " << r(j);
//                }
//                for(k=0;k<dim;k++){
//                    myfile << " " << v(k);
//                }

//                myfile << "\n";
//            }

            v_half = v + F*dt/2;

            if(can_move==true && removed_atom==false){
                r_new = r + v_half*dt;
            } else {
                r_new = r;
            }


            if(t==202 && can_move==true && removed_atom==false){  // Counting atoms
                moving_atoms = moving_atoms + 1;
            }
            if(t==202 && can_move==false){
                non_moving_atoms = non_moving_atoms + 1;
            }

            if(t==201 && can_move==true){   // Making a fluid with half the density
                tall = randu();
                if(tall<0.5){
                    atoms[i].remove(true);
                    removed_atoms = removed_atoms + 1;
                }
            }

            for(j=0;j<dim;j++){  // periodic boundary conditions

                r_new(j) = fmod(r_new(j), (L));

                if (r_new(j)<0){
                    r_new(j) = L + r_new(j);
                }
            }

            if(t>(n-1000) && r(2)<L/2 && r_new(2)>L/2 && v_half(2)>0){   // Finding flow through z = L/2
                U = U + 1.0;
            }
            if(t>(n-1000) && r(2)>L/2 && r_new(2)<L/2 && v_half(2)<0){
                U = U - 1.0;
            }

            if(t>(n-1000) && can_move==true && removed_atom==false){    // Finding u(r)
                r_xy = sqrt((r(0)-L/2)*(r(0)-L/2) + (r(1)-L/2)*(r(1)-L/2));
                v_xyz = sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2));

                for(j=0;j<n_u;j++){
                    if(radius(j) < r_xy && r_xy < radius(j+1)){
                        u(j) = u(j) + v_xyz;
                        uz(j) = uz(j) + v(2);
                        number_of_elements(j) = number_of_elements(j) + 1.0;
                    }
                }
            }

            atoms[i].setPosition(r_new);

        }


        if(t>(n-1000)){
            for(i=0;i<n_u;i++){
                u(i) = u(i)/number_of_elements(i);                
                myfile << u(i) << " ";
                u(i) = 0.0;
            }
            for(i=0;i<n_u;i++){
                uz(i) = uz(i)/number_of_elements(i);
                myfile << uz(i) << " ";
                uz(i) = 0.0;
            }
            for(i=0;i<n_u;i++){
                density(i) = number_of_elements(i)/volume_cyl(i);
                myfile << density(i) << " ";
                number_of_elements(i) = 0.0;
            }

            myfile << U << "\n"; // Finding flow
            U = 0;
        }

        if(t==200){
            Cylinder Cyl;

            for(i=0;i<numpart;i++){
                Cyl.testCylinder(&atoms[i],R,L);
            }
        }

//        if(t==200){
//            Spheres Sph(spheres,L);
//            for(i=0;i<numpart;i++){
//                Sph.findSpheres(&atoms[i],sphere_pos_x,sphere_pos_y,sphere_pos_z,sphere_radius);
//                if(i==0){
//                    sphere_volume = Sph.findSphereVolume(pi,sphere_radius);
//                }
//            }
//        }

//        E_tot = E_k + pot;

//        T = 2*E_k/(3.0*numpart);
//        P = numpart/V*T + 1/(3*V)*Fr_sum;

//        if(t>201){
//            T = 2*E_k/(3.0*moving_atoms);
//            P = moving_atoms/(pi*R*R*L)*T + 1.0/(3*pi*R*R*L)*Fr_sum;
//            P = moving_atoms/(V-sphere_volume)*T + 1/(3*(V-sphere_volume))*Fr_sum;
//            myfile << T << " " << P << "\n";
//        }

//        if(t>300){
//            mean_disp = 1.0/moving_atoms*r_sum;
//            myfile << mean_disp << "\n";
//        }



        //        for(i=0;i<numpart;i++){  // Berendsen thermostat
        //            v = atoms[i].getVelocity();
        //            gamma = sqrt(1 + dt/tau*(T_bath/T - 1));
        //            for(j=0;j<dim;j++){
        //                v(j) = v(j)*gamma;
        //            }
        //            atoms[i].setVelocity(v);
        //        }


        if(t<170){
            for(i=0;i<numpart;i++){  // Andersen thermostat
                tall = randu();
                if(tall<(dt/tau)){
                    v_new = randn(3)*sqrt(T_thermalize); //*sqrt(T_bath);
                    atoms[i].setVelocity(v_new);
                }
            }
        }

        if(t>205 && t<500){
            for(i=0;i<numpart;i++){  // Andersen thermostat
                can_move = atoms[i].getMove();
                removed_atom = atoms[i].testRemoved();
                tall = randu();
                if(tall<(dt/tau) && can_move==true && removed_atom==false){
                    v_new = randn(3)*sqrt(T_thermalize2); //*sqrt(T_bath);
                    atoms[i].setVelocity(v_new);
                }
            }
        }


        E_k = 0.0;
        pot = 0.0;
        Fr_sum = 0.0;
        r_sum = 0.0;

        Cell_container Con(cells_x);


        for(i=0;i<numpart;i++){  // Setter nye posisjoner i celler
            r = atoms[i].getPosition();
            nx = int (r(0)/r_cut);
            ny = int (r(1)/r_cut);
            nz = int (r(2)/r_cut);

            if(nx>cells_x-1){
                nx = cells_x-1;
            }
            if(ny>cells_x-1){
                ny = cells_x-1;
            }
            if(nz>cells_x-1){
                nz = cells_x-1;
            }

            Con.container(nx,ny,nz,&atoms[i]);

        }



        for(i=0;i<cells_x;i++){  // Finn kraft og fart
            for(j=0;j<cells_x;j++){
                for(k=0;k<cells_x;k++){
                    celle = Con.getContainer(i,j,k);

                    for(l=0;l<celle.size();l++){

                        removed_atom = celle[l]->testRemoved();
                        can_move = celle[l]->getMove();

                        if(removed_atom==false && can_move==true){

                            r = celle[l]->getPosition();

                            for(g=0;g<dim;g++){
                                F(g) = 0;
                            }

                            for(g=0;g<celle.size();g++){   // egen celle
                                if(g!=l){
                                    removed_atom2 = celle[g]->testRemoved();
                                    can_move2 = celle[g]->getMove();

                                    if(removed_atom2==false){
                                        r2 = celle[g]->getPosition();

                                        d = r - r2;
                                        d_length = sqrt(d(0)*d(0) + d(1)*d(1) + d(2)*d(2));
                                        pot = pot + 4*(pow (1/d_length,12) - pow (1/d_length,6));

                                        for(h=0;h<dim;h++){
                                            F_bidrag(h) = 24*(2*pow (1/d_length,12) - pow (1/d_length,6))*d(h)/(d_length*d_length);
                                            F(h) = F(h) + F_bidrag(h);
                                            if(can_move2==true){
                                                Fr_sum = Fr_sum + 0.5*F_bidrag(h)*d(h);
                                            } else{
                                                Fr_sum = Fr_sum + F_bidrag(h)*d(h);
                                            }
                                        }
                                    }
                                }
                            }


                            for(g=-1;g<2;g++){   // naboceller
                                for(h=-1;h<2;h++){
                                    for(o=-1;o<2;o++){
                                        if(!(g==0 && h==0 && o==0)){

                                            i2 = i+g;
                                            j2 = j+h;
                                            k2 = k+o;

                                            if((i+g) < 0){
                                                i2 = cells_x-1;
                                            }
                                            if((j+h) < 0){
                                                j2 = cells_x-1;
                                            }
                                            if((k+o) < 0){
                                                k2 = cells_x-1;
                                            }

                                            if((i+g) == cells_x){
                                                i2 = 0;
                                            }
                                            if((j+h) == cells_x){
                                                j2 = 0;
                                            }
                                            if((k+o) == cells_x){
                                                k2 = 0;
                                            }

                                            celle2 = Con.getContainer(i2,j2,k2);


                                            for(p=0;p<celle2.size();p++){

                                                removed_atom2 = celle2[p]->testRemoved();
                                                can_move2 = celle2[p]->getMove();

                                                if(removed_atom2==false){

                                                    r2 = celle2[p]->getPosition();
                                                    d = r-r2;

                                                    for(q=0;q<dim;q++){    // minimum image convention
                                                        if(fabs(d(q)-L) < fabs(d(q))){
                                                            if(fabs(d(q)-L) < fabs(d(q)+L)){
                                                                r_min(q) = d(q)-L;
                                                            }
                                                        }
                                                        if(fabs(d(q)) < fabs(d(q)+L)){
                                                            if(fabs(d(q)) < fabs(d(q)-L)){
                                                                r_min(q) = d(q);
                                                            }
                                                        }
                                                        if(fabs(d(q)+L) < fabs(d(q)-L)){
                                                            if(fabs(d(q)+L) < fabs(d(q))){
                                                                r_min(q) = d(q)+L;
                                                            }
                                                        }
                                                    }


                                                    d_length = sqrt(r_min(0)*r_min(0) + r_min(1)*r_min(1) + r_min(2)*r_min(2));

                                                    pot = pot + 4*(pow (1/d_length,12) - pow (1/d_length,6));

                                                    for(q=0;q<dim;q++){
                                                        F_bidrag(q) = 24*(2*pow (1/d_length,12) - pow (1/d_length,6))*r_min(q)/(d_length*d_length);
                                                        F(q) = F(q) + F_bidrag(q);
                                                        if(can_move2==true){
                                                            Fr_sum = Fr_sum + 0.5*F_bidrag(q)*r_min(q);
                                                        } else {
                                                            Fr_sum = Fr_sum + F_bidrag(q)*r_min(q);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }


                            if(t>300){ // && t<1300){
                                F(2) = F(2) + 0.1; // 0.1;  //F_x
                            }

                            v = celle[l]->getVelocity();
                            F_old = celle[l]->getForce();

                            v_half = v + F_old*dt/2;
                            v_new = v_half + F*dt/2;                            

                            celle[l]->setForce(F);
                            celle[l]->setVelocity(v_new);

//                            E_k = E_k + 1.0/2.0*(v_new(0)*v_new(0) + v_new(1)*v_new(1) + v_new(2)*v_new(2));
//                            moving_atoms2 = moving_atoms2 + 1;

                        }
                    }
//                    if(t>300){     // Calculating pressure in each cell

//                        if(i==(cells_x-1) && j==(cells_x-1) && k==(cells_x-1)){
//                            V_cell = pow (L - (cells_x-1)*r_cut,3);
//                        }
//                        if((i==(cells_x-1) && j==(cells_x-1) && k!=(cells_x-1)) || (i==(cells_x-1) && j!=(cells_x-1) && k==(cells_x-1)) || (i!=(cells_x-1) && j==(cells_x-1) && k==(cells_x-1))){
//                            V_cell = pow (L - (cells_x-1)*r_cut,2)*r_cut;
//                        }
//                        if((i==(cells_x-1) && j!=(cells_x-1) && k!=(cells_x-1)) || (i!=(cells_x-1) && j!=(cells_x-1) && k==(cells_x-1)) || (i!=(cells_x-1) && j==(cells_x-1) && k!=(cells_x-1))){
//                            V_cell = (L - (cells_x-1)*r_cut)*pow (r_cut,2);
//                        }
//                        else{
//                            V_cell = pow (r_cut,3);
//                        }

//                        if(moving_atoms2 == 0){
//                            T = 0;
//                        }
//                        else{
//                            T = 2.0*E_k/(3.0*moving_atoms2);
//                        }
//                        P = (float (moving_atoms2))/V_cell*T + 1.0/(3.0*V_cell)*Fr_sum;

////                        myfile << P << " ";
//                    }

//                    E_k = 0;
//                    Fr_sum = 0;
//                    moving_atoms2 = 0;
                }
            }
        }

//        if(t>300){
//            myfile << "\n";
//        }
    }
    myfile.close();

//    double porosity = (float (moving_atoms))/numpart;
//    cout << "porosity = " << porosity << endl;

//    //sphere_volume = Sph.findSphereVolume(pi);
//    double porosity2 = (V - sphere_volume)/V;
//    cout << "theoretical porosity = " << porosity2 << endl;
}

//porosity = 0.404344
//theoretical porosity = 0.0469301

