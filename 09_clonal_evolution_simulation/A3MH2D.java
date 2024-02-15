package Model;

import Framework.GridsAndAgents.AgentGrid2D;
import Framework.GridsAndAgents.AgentSQ2Dunstackable;
import Framework.Gui.*;
import Framework.Tools.FileIO;
import Framework.Rand;
import org.lwjgl.Sys;

import static Framework.Util.*;
import java.util.ArrayList;

class A3Parameters2D {
    public double mu = 1e-8;
    public double mua3 = 0;
    public double delta = 0;
    public double stg = 0;

    public double sn = 1e-3;
    public double sp = 0.1;

    public double Tm = 5e6;

    public double a3l = 0;
    public double a3m = 1;
    public double a3h = 0;

    public double pa3l = 0.001;
    public double pa3m = 0.01;
    public double pa3h = 0.1;

    public double Ta3 = (int)(Tm*a3l*pa3l + Tm*a3m*pa3m + Tm*a3h*pa3h);

    public double pNeg = 0.4;
    public double pNeu = 0.5;
    public double pPos = 0.1;

    public double a3lpNegna3 = 0.1;
    public double a3lpNeuna3 = 0.5;
    public double a3lpPosna3 = 0.4;
    public double a3lpdelta = 0;
    public double a3lpNega3 = 0.1;
    public double a3lpNeua3 = 0.5;
    public double a3lpPosa3 = 0.4;

    public double a3mpNegna3 = pNeg;
    public double a3mpNeuna3 = pNeu;
    public double a3mpPosna3 = pPos;
    public double a3mpdelta = 0;
    public double a3mpNega3 = pNeg;
    public double a3mpNeua3 = pNeu;
    public double a3mpPosa3 = pPos;

    public double a3hpNegna3 = 0.1;
    public double a3hpNeuna3 = 0.5;
    public double a3hpPosna3 = 0.4;
    public double a3hpdelta = 0;
    public double a3hpNega3 = 0.1;
    public double a3hpNeua3 = 0.5;
    public double a3hpPosa3 = 0.4;

    public double birth_rate = 0.5;
    public double death_rate = 0.5;
    public int r0 = 10;
    int sideLen = 200;

    // tracking variables
    public int KnMAX = 0;
    public int Knna3MAX = 0;
    public int Kna3MAX = 0;

    public int KoMAX = 0;
    public int Kona3MAX = 0;
    public int Koa3MAX = 0;

    public int KpMAX = 0;
    public int Kpna3MAX = 0;
    public int Kpa3MAX = 0;

    int delete_thresh = 75; // ignore clone sizes smaller than this in Muller plots


}

class A3Cell2D extends AgentSQ2Dunstackable<A3MH2D> {
    int kn;
    int ko;
    int kp;

    int kna3;
    int koa3;
    int kpa3;

    int knna3;
    int kona3;
    int kpna3;

    double fit;

    int progenyID;
    int parentID;



    A3Cell2D Init(int kn0, int ko0, int kp0, int kna30, int koa30, int kpa30,
                  int knna30, int kona30, int kpna30, int progenyID0, int parentID0){
        kn = kn0;
        ko = ko0;
        kp = kp0;
        kna3 = kna30;
        koa3 = koa30;
        kpa3 = kpa30;
        knna3 = knna30;
        kona3 = kona30;
        kpna3 = kpna30;
        progenyID = progenyID0;
        parentID = parentID0;
        fit = (Math.pow(1.0+G.sp,(double)kp)/Math.pow(1.0+G.sn,(double)kn));
        return this;
    }

    A3Cell2D Mutate(){
        boolean mutated = false;

        boolean neg_mutated = false;
        boolean neu_mutated = false;
        boolean pos_mutated = false;

        boolean na3_neg_mutated = false;
        boolean na3_neu_mutated = false;
        boolean na3_pos_mutated = false;

        boolean a3_neg_mutated = false;
        boolean a3_neu_mutated = false;
        boolean a3_pos_mutated = false;

        double phita3;

        double pa3Negna3;
        double pa3Neuna3;
        double pa3Posna3;

        double pa3Nega3;
        double pa3Neua3;
        double pa3Posa3;

        // random mutation
        if (G.rn.Double() < ( G.Tm * G.mu)) {
            mutated = true;
            double n = G.rn.Double();
            if (n < G.a3l) {
                // if mutation hits low A3 motifs gene
                phita3 = G.pa3l;

                pa3Negna3 = G.a3lpNegna3;
                pa3Neuna3 = G.a3lpNeuna3;
                pa3Posna3 = G.a3lpPosna3;

                pa3Nega3 = G.a3lpNega3;
                pa3Neua3 = G.a3lpNeua3;
                pa3Posa3 = G.a3lpPosa3;

            } else if (n < G.a3l + G.a3m) {
                // if mutation hits moderate A3 motifs gene
                phita3 = G.pa3m;

                pa3Negna3 = G.a3mpNegna3;
                pa3Neuna3 = G.a3mpNeuna3;
                pa3Posna3 = G.a3mpPosna3;

                pa3Nega3 = G.a3mpNega3;
                pa3Neua3 = G.a3mpNeua3;
                pa3Posa3 = G.a3mpPosa3;
            } else {
                // if mutation hits high A3 motifs gene
                phita3 = G.pa3h;

                pa3Negna3 = G.a3hpNegna3;
                pa3Neuna3 = G.a3hpNeuna3;
                pa3Posna3 = G.a3hpPosna3;

                pa3Nega3 = G.a3hpNega3;
                pa3Neua3 = G.a3hpNeua3;
                pa3Posa3 = G.a3hpPosa3;
            }

            if (G.rn.Double() < phita3) {
                // if mutation hits A3 motif
                n = G.rn.Double();
                if (n < pa3Nega3) {
                    kn++;
                    kna3++;
                    if (kn > G.KnMAX) { G.KnMAX++; }
                    if (kna3 > G.Kna3MAX) { G.Kna3MAX++; }
                } else if (n < pa3Nega3 + pa3Neua3) {
                    if (koa3 > 0) {
                        if (G.rn.Double() < ko*G.stg) {
                            if (ko == G.KoMAX) {
                                int cntKoMAX = 0;
                                for (int i = 0; i < G.sideLen*G.sideLen; i++) {
                                    if (G.GetAgent(i) != null) {
                                        if (G.GetAgent(i).ko == G.KoMAX) {
                                            cntKoMAX++;
                                        }

                                    }
                                }

                                if (cntKoMAX == 1) {G.KoMAX--;}

                            } else if (koa3 == G.Koa3MAX) {
                                int cntKoa3MAX = 0;
                                for (int i = 0; i < G.sideLen*G.sideLen; i++) {
                                    if (G.GetAgent(i) != null) {

                                        if (G.GetAgent(i).koa3 == G.Koa3MAX) {
                                            cntKoa3MAX++;
                                        }
                                    }
                                }

                                if (cntKoa3MAX == 1) {G.Koa3MAX--;}
                            }
                            ko--;
                            koa3--;
                            kp++;
                            kpa3++;
                            if (kp > G.KpMAX) { G.KpMAX++; }
                            if (kpa3 > G.Kpa3MAX) { G.Kpa3MAX++; }

                        } else {
                            ko++;
                            koa3++;

                            if (ko > G.KoMAX) { G.KoMAX++; }
                            if (koa3 > G.Koa3MAX) { G.Koa3MAX++; }
                        }
                    } else {
                        ko++;
                        koa3++;

                        if (ko > G.KoMAX) { G.KoMAX++; }
                        if (koa3 > G.Koa3MAX) { G.Koa3MAX++; }
                    }
                } else {
                    kp++;
                    kpa3++;
                    if (kp > G.KpMAX) { G.KpMAX++; }
                    if (kpa3 > G.Kpa3MAX) { G.Kpa3MAX++; }
                }

            } else {
                // if mutation hits non A3 motif
                n = G.rn.Double();
                if (n < pa3Negna3) {
                    kn++;
                    knna3++;
                    if (kn > G.KnMAX) { G.KnMAX++; }
                    if (knna3 > G.Knna3MAX) { G.Knna3MAX++; }
                } else if (n < pa3Negna3 + pa3Neuna3) {
                    ko++;
                    kona3++;
                    if (ko > G.KoMAX) { G.KoMAX++; }
                    if (kona3 > G.Kona3MAX) { G.Kona3MAX++; }
                } else {
                    kp++;
                    kpna3++;
                    if (kp > G.KpMAX) { G.KpMAX++; }
                    if (kpna3 > G.Kpna3MAX) { G.Kpna3MAX++; }
                }
            }
        }

        // additional a3 mutation
        if (G.mua3 != 0) {
            if (G.rn.Double() < ( G.Ta3 * G.mua3)) {
                mutated = true;
                double n = G.rn.Double();
                if (n < G.a3l) {
                    // if mutation hits low A3 motifs gene
                    pa3Nega3 = G.a3lpNega3;
                    pa3Neua3 = G.a3lpNeua3;
                    pa3Posa3 = G.a3lpPosa3;

                } else if (n < G.a3l + G.a3m) {
                    // if mutation hits moderate A3 motifs gene
                    pa3Nega3 = G.a3mpNega3;
                    pa3Neua3 = G.a3mpNeua3;
                    pa3Posa3 = G.a3mpPosa3;
                } else {
                    // if mutation hits high A3 motifs gene
                    pa3Nega3 = G.a3hpNega3;
                    pa3Neua3 = G.a3hpNeua3;
                    pa3Posa3 = G.a3hpPosa3;
                }

                n = G.rn.Double();
                if (n < pa3Nega3) {
                    kn++;
                    kna3++;
                    if (kn > G.KnMAX) { G.KnMAX++; }
                    if (kna3 > G.Kna3MAX) { G.Kna3MAX++; }
                } else if (n < pa3Nega3 + pa3Neua3) {
                    if (koa3 > 0) {
                        if (G.rn.Double() < ko*G.stg) {
                            if (ko == G.KoMAX) {
                                int cntKoMAX = 0;
                                for (int i = 0; i < G.sideLen*G.sideLen; i++) {
                                    if (G.GetAgent(i) != null) {
                                        if (G.GetAgent(i).ko == G.KoMAX) {
                                            cntKoMAX++;
                                        }

                                    }
                                }

                                if (cntKoMAX == 1) {G.KoMAX--;}

                            } else if (koa3 == G.Koa3MAX) {
                                int cntKoa3MAX = 0;
                                for (int i = 0; i < G.sideLen*G.sideLen; i++) {
                                    if (G.GetAgent(i) != null) {

                                        if (G.GetAgent(i).koa3 == G.Koa3MAX) {
                                            cntKoa3MAX++;
                                        }
                                    }
                                }

                                if (cntKoa3MAX == 1) {G.Koa3MAX--;}
                            }
                            ko--;
                            koa3--;
                            kp++;
                            kpa3++;
                            if (kp > G.KpMAX) { G.KpMAX++; }
                            if (kpa3 > G.Kpa3MAX) { G.Kpa3MAX++; }

                        } else {
                            ko++;
                            koa3++;

                            if (ko > G.KoMAX) { G.KoMAX++; }
                            if (koa3 > G.Koa3MAX) { G.Koa3MAX++; }
                        }
                    } else {
                        ko++;
                        koa3++;

                        if (ko > G.KoMAX) { G.KoMAX++; }
                        if (koa3 > G.Koa3MAX) { G.Koa3MAX++; }
                    }


                } else {
                    kp++;
                    kpa3++;
                    if (kp > G.KpMAX) { G.KpMAX++; }
                    if (kpa3 > G.Kpa3MAX) { G.Kpa3MAX++; }
                }
            }
        }

        if (mutated) {
            parentID = progenyID;
            progenyID = G.progenyNextID;
            G.progenyToParentIDs[progenyID] = parentID;

            G.all_status[progenyID] = kn + ko + kp;

            G.neg_status[progenyID] = kn;
            G.neu_status[progenyID] = ko;
            G.pos_status[progenyID] = kp;

            G.nega3_status[progenyID] = kna3;
            G.neua3_status[progenyID] = koa3;
            G.posa3_status[progenyID] = kpa3;

            G.negna3_status[progenyID] = knna3;
            G.neuna3_status[progenyID] = kona3;
            G.posna3_status[progenyID] = kpna3;

            G.progenyNextID++;

            fit = (Math.pow(1.0+G.sp,(double)kp)/Math.pow(1.0+G.sn,(double)kn));
            G.genotype_status[progenyID] = kp*(G.KnMAX+1)*(G.KoMAX+1) + kn*(G.KoMAX+1) + ko;
            G.fitness_status[progenyID] = kp*(G.KnMAX+1) + kn;
        }

        return this;
    }

    A3Cell2D Divide(){
        int nDivOptions = G.MapEmptyHood(G.neighborhood,Xsq(),Ysq());
        if(nDivOptions==0){
            return null;
        }
        int nextAgentID = G.neighborhood[G.rn.Int(nDivOptions)];

        double first = (nextAgentID/G.sideLen) - (G.sideLen/2);
        double second = (nextAgentID%G.sideLen) - (G.sideLen/2);

        if ((first*first + second*second) <= (G.confRadius*G.confRadius)) {
            return G.NewAgentSQ(nextAgentID).Init(this.kn,this.ko,this.kp,this.kna3,this.koa3,this.kpa3,this.knna3,this.kona3,this.kpna3,this.progenyID,this.parentID).Mutate();
        } else {
            return this;
        }

    }

    void Step(){
        // Passengers lower birth rate; Drivers raise birth rate
        if(G.rn.Double()<(Math.pow(1.0+G.sp,(double)kp)/Math.pow(1.0+G.sn,(double)kn)*G.birth_rate)){
            Divide();
        }

        // constant death rate
        if(G.rn.Double()<(G.death_rate )){
            Dispose();
            return;
        }

    }
}

public class A3MH2D extends AgentGrid2D<A3Cell2D> {

    // arrays used to store parentIDs for building Muller plots (and associated driver numbers)
    public int[] progenyToParentIDs = new int[7000000];

    public int [] all_status = new int[7000000];

    public int[] neg_status = new int[7000000];
    public int[] neu_status = new int[7000000];
    public int[] pos_status = new int[7000000];

    public int[] nega3_status = new int[7000000];
    public int[] neua3_status = new int[7000000];
    public int[] posa3_status = new int[7000000];

    public int[] negna3_status = new int[7000000];
    public int[] neuna3_status = new int[7000000];
    public int[] posna3_status = new int[7000000];

    public int[] genotype_status = new int[7000000];
    public int[] fitness_status = new int[7000000];

    // parameters

    public double mu;
    public double mua3;
    public double delta;
    public double stg;

    public double sn;
    public double sp;

    public double Tm;

    public double a3l;
    public double a3m;
    public double a3h;

    public double pa3l;
    public double pa3m;
    public double pa3h;

    public double pNeg;
    public double pNeu;
    public double pPos;

    public double Ta3;

    public double a3lpNegna3;
    public double a3lpNeuna3;
    public double a3lpPosna3;
    public double a3lpdelta;
    public double a3lpNega3;
    public double a3lpNeua3;
    public double a3lpPosa3;

    public double a3mpNegna3;
    public double a3mpNeuna3;
    public double a3mpPosna3;
    public double a3mpdelta;
    public double a3mpNega3;
    public double a3mpNeua3;
    public double a3mpPosa3;

    public double a3hpNegna3;
    public double a3hpNeuna3;
    public double a3hpPosna3;
    public double a3hpdelta;
    public double a3hpNega3;
    public double a3hpNeua3;
    public double a3hpPosa3;

    public int confRadius;
    public double birth_rate;
    public double death_rate;
    public int r0;
    int sideLen;

    // tracking variables
    public int KnMAX;
    public int Knna3MAX;
    public int Kna3MAX;

    public int KoMAX;
    public int Kona3MAX;
    public int Koa3MAX;

    public int KpMAX;
    public int Kpna3MAX;
    public int Kpa3MAX;

    public int progenyNextID = 2; // assumes model is initialized with all cell's progenyID = 1

    // neighborhoods
    int[]neighborhood = MooreHood(false);

    Rand rn0 = new Rand();
    public int seed = rn0.Int(10000);
    Rand rn = new Rand(seed);

    // constructor from parameters
    A3MH2D(A3Parameters2D p){
        super(p.sideLen,p.sideLen, A3Cell2D.class,false,false);
        for (int x = 0; x < p.r0; x++) {
            for (int y = 0; y < p.r0; y++) {
                if ((x-p.r0/2)*(x-p.r0/2)+(y-p.r0/2)*(y-p.r0/2) < p.r0*p.r0/4) {
                    NewAgentSQ(x + p.sideLen / 2 - p.r0 / 2, y + p.sideLen / 2 - p.r0 / 2).Init(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0);
                }
            }
        }

        this.mu = p.mu;
        this.mua3 = p.mua3;
        this.delta = p.delta;

        this.stg = p.stg;

        this.sn = p.sn;
        this.sp = p.sp;

        this.Tm = p.Tm;

        this.a3l = p.a3l;
        this.a3m = p.a3m;
        this.a3h = p.a3h;

        this.Ta3 = p.Ta3;

        this.pa3l = p.pa3l;
        this.pa3m = p.pa3m;
        this.pa3h = p.pa3h;

        this.pNeg = p.pNeg;
        this.pNeu = p.pNeu;
        this.pPos = p.pPos;

        this.a3lpNegna3 = p.a3lpNegna3;
        this.a3lpNeuna3 = p.a3lpNeuna3;
        this.a3lpPosna3 = p.a3lpPosna3;
        this.a3lpdelta = p.a3lpdelta;
        this.a3lpNega3 = p.a3lpNega3;
        this.a3lpNeua3 = p.a3lpNeua3;
        this.a3lpPosa3 = p.a3lpPosa3;

        this.a3mpNegna3 = p.a3mpNegna3;
        this.a3mpNeuna3 = p.a3mpNeuna3;
        this.a3mpPosna3 = p.a3mpPosna3;
        this.a3mpdelta = p.a3mpdelta;
        this.a3mpNega3 = p.a3mpNega3;
        this.a3mpNeua3 = p.a3mpNeua3;
        this.a3mpPosa3 = p.a3mpPosa3;

        this.a3hpNegna3 = p.a3hpNegna3;
        this.a3hpNeuna3 = p.a3hpNeuna3;
        this.a3hpPosna3 = p.a3hpPosna3;
        this.a3hpdelta = p.a3hpdelta;
        this.a3hpNega3 = p.a3hpNega3;
        this.a3hpNeua3 = p.a3hpNeua3;
        this.a3hpPosa3 = p.a3hpPosa3;

        this.sideLen = p.sideLen;
        this.r0 = p.r0;

        this.birth_rate = p.birth_rate;
        this.death_rate = p.death_rate;
        this.confRadius = p.sideLen/2;

        // iterators / indicators
        this.KnMAX = p.KnMAX;
        this.Knna3MAX = p.Knna3MAX;
        this.Kna3MAX = p.Kna3MAX;

        this.KoMAX = p.KoMAX;
        this.Kona3MAX = p.Kona3MAX;
        this.Koa3MAX = p.Koa3MAX;

        this.KpMAX = p.KpMAX;
        this.Kpna3MAX = p.Kpna3MAX;
        this.Kpa3MAX = p.Kpa3MAX;
    }

    // Step function for PD model ("steps" all cells through birth/death/mutation)
    void OriginalStep(){
        for (A3Cell2D c:this) {
            c.Step();
        }
        CleanShuffle(rn);
    }

    public static int SingleSim(A3Parameters2D p, int totalTime, int modifier, int sim, boolean saveTimelinesBoolean, String foldername) {
        A3MH2D model = new A3MH2D(p);

        // VISUALIZE
        UIWindow win = new UIWindow("Window seed: " + model.seed, true);
        UIGrid Vis = new UIGrid(model.sideLen, model.sideLen, 3);
        TickTimer tt = new TickTimer();
        win.AddCol(0, Vis);
        win.RunGui();

        // INITIALIZE ARRAY
        // 0. time 1. pop
        // 2-10. neg, neu, pos, neg a3, neu a3, pos a3, neg na3, neu na3, pos na3 diversity
        // 11-19. neg, neu, pos, neg a3, neu a3, pos a3, neg na3, neu na3, pos na3 n
        // 22-28. neg, neu, pos, neg a3, neu a3, pos a3, neg na3, neu na3, pos na3 max
        // 29. fitness diversity
        // 30. genotype diversity
        // 31. number of genotypes
        // fill in time
        double[][] everything = new double[32][(totalTime/modifier)+1];
        for(int iii=0; iii<32; iii++) {
            for(int jjj=0; jjj<(totalTime/modifier)+1; jjj++) {

                if (iii == 0) {
                    everything[iii][jjj] = (double)modifier * jjj;
                } else {
                    everything[iii][jjj] = 0.0;
                }
            }
        }

        int expectedProgeny = model.progenyToParentIDs.length;
        int[][] mullerGenetic = new int[expectedProgeny][(totalTime/modifier)+1]; // track the first [expectedProgeny] genetic clones
        for (int jjj = 0; jjj < (totalTime / modifier)+1; jjj++) {
            for (int iii = 0; iii < expectedProgeny; iii++) {
                mullerGenetic[iii][jjj] = 0;
            }
        }

        // new gif
        String baseFilename =  foldername +
                "br" + (int) (model.birth_rate*100) +
                "dr" + (int) (model.death_rate*100) +
                "a3l" + (int) (model.a3l*100 ) +
                "a3m" + (int) (model.a3m*100 ) +
                "a3h" + (int) (model.a3h*100 ) +
                "seed" + model.seed +
                "sim" + sim;
        FileIO everythingOutputFile =  new FileIO((baseFilename + ".csv"),( (sim == 0) ? "w" : "a" ));
        String gifFilename = baseFilename + ".gif";
        GifMaker myGif = new GifMaker(gifFilename, 1,true);

        int j = 0;
        boolean keepStepping = true;
        for (int i = 0; i < totalTime+1; i++) {
            // episodic APOBEC3
            if (i == 0){
                model.mua3 = model.mua3 + model.delta;
            }
            if (i == 50){
                model.mua3 = model.mua3 - model.delta;
            }

            DrawCells(model, Vis);
            tt.TickPause(1);
            myGif.AddFrame(Vis);
            System.out.println("time: " + i + " pop: " + model.Pop() + " max mut: " + (model.KnMAX + model.KoMAX + model.KpMAX));

            if (i != 0) {
                if (keepStepping) { model.OriginalStep(); }
            }

            if (model.Pop() == 0) { keepStepping = false; }

            // WRITE OUT VALUES
            if (i % modifier == 0) {

                if (saveTimelinesBoolean) {
                    everything[1][j] = model.Pop(); // total pop size

                    everything[2][j] = GetDiversity(model, 0); // neg div
                    everything[3][j] = GetDiversity(model, 1); // neu div
                    everything[4][j] = GetDiversity(model, 2); // pos div

                    everything[5][j] = GetDiversity(model, 3); // neg a3 div
                    everything[6][j] = GetDiversity(model, 4); // neu a3 div
                    everything[7][j] = GetDiversity(model, 5); // pos a3 div

                    everything[8][j] = GetDiversity(model, 6); // neg na3 div
                    everything[9][j] = GetDiversity(model, 7); // neu na3 div
                    everything[10][j] = GetDiversity(model, 8); // pos na3 div

                    everything[11][j] = CountMutators(model, 0); // neg num
                    everything[12][j] = CountMutators(model, 1); // neu num
                    everything[13][j] = CountMutators(model, 2); // pos num

                    everything[14][j] = CountMutators(model, 3); // neg a3 num
                    everything[15][j] = CountMutators(model, 4); // neu a3 num
                    everything[16][j] = CountMutators(model, 5); // pos a3 num

                    everything[17][j] = CountMutators(model, 6); // neg na3 num
                    everything[18][j] = CountMutators(model, 7); // neu na3 num
                    everything[19][j] = CountMutators(model, 8); // pos na3 num

                    everything[20][j] = model.KnMAX; // neg max
                    everything[21][j] = model.KoMAX; // neu max
                    everything[22][j] = model.KpMAX; // pos max

                    everything[23][j] = model.Kna3MAX; // neg a3 max
                    everything[24][j] = model.Koa3MAX; // neu a3 max
                    everything[25][j] = model.Kpa3MAX; // pos a3 max

                    everything[26][j] = model.Knna3MAX; // neg na3 max
                    everything[27][j] = model.Kona3MAX; // neu na3 max
                    everything[28][j] = model.Kpna3MAX; // pos na3 max

                    everything[29][j] = GetFitDiversity(model);
                    everything[30][j] = GetGenotypeDiversity(model);
                    everything[31][j] = GetNumGenotypes(model);

                    // add all the Muller information
                    for (int k = 0; k < (model.sideLen*model.sideLen); k++) {
                        A3Cell2D c = model.GetAgent(k);
                        if (c != null) {
                            mullerGenetic[c.progenyID][j] = mullerGenetic[c.progenyID][j] + 1;
                        }
                    }
                }

                j++;
            }
        }

        if (saveTimelinesBoolean) {

            System.out.println("Building everything to output...");

            // OUTPUT "EVERYTHING" ARRAYS TO FILE
            BuildEverythingSingle(everythingOutputFile, everything, totalTime, modifier);

            FileIO parentsFile = new FileIO((baseFilename + "_parents.csv"), "w");
            FileIO geneticMullerFile = new FileIO((baseFilename + "_clones.csv"), "w");
            FileIO allStatusFile = new FileIO((baseFilename + "_allStatus.csv"), "w");

            FileIO negStatusFile = new FileIO((baseFilename + "_negStatus.csv"), "w");
            FileIO neuStatusFile = new FileIO((baseFilename + "_neuStatus.csv"), "w");
            FileIO posStatusFile = new FileIO((baseFilename + "_posStatus.csv"), "w");

            FileIO nega3StatusFile = new FileIO((baseFilename + "_nega3Status.csv"), "w");
            FileIO neua3StatusFile = new FileIO((baseFilename + "_neua3Status.csv"), "w");
            FileIO posa3StatusFile = new FileIO((baseFilename + "_posa3Status.csv"), "w");

            FileIO negna3StatusFile = new FileIO((baseFilename + "_negna3Status.csv"), "w");
            FileIO neuna3StatusFile = new FileIO((baseFilename + "_neuna3Status.csv"), "w");
            FileIO posna3StatusFile = new FileIO((baseFilename + "_posna3Status.csv"), "w");

            FileIO genotypeStatusFile = new FileIO((baseFilename + "_genotypeStatus.csv"), "w");
            FileIO fitnessStatusFile = new FileIO((baseFilename + "_fitnessStatus.csv"), "w");


            SaveAll(parentsFile, geneticMullerFile, allStatusFile,
                    negStatusFile, neuStatusFile, posStatusFile,
                    nega3StatusFile, neua3StatusFile, posa3StatusFile,
                    negna3StatusFile, neuna3StatusFile, posna3StatusFile,
                    genotypeStatusFile, fitnessStatusFile,
                    mullerGenetic, model.progenyToParentIDs, model.all_status,
                    model.neg_status, model.neu_status, model.pos_status,
                    model.nega3_status, model.neua3_status, model.posa3_status,
                    model.negna3_status, model.neuna3_status, model.posna3_status,
                    model.genotype_status, model.fitness_status,
                    expectedProgeny, model.progenyNextID, totalTime, modifier, p.delete_thresh);

            System.out.println("Simulation finished...");
            System.out.println("Okay to close...");

        }

        //System.out.println(model.progenyNextID);

        //win.Close();

        return model.seed;
    }

    public static void main(String[] args) {
        /*
            SINGLE SIMULATION, constrained to a circular domain
                - uncomment out the following section to run a single sim
                - outputs basic information ("everything" array) and GIFs and Muller plot information
         */

        int modifier = 10; // when to save data
        int totalTime = 600; // total simulation time
        int nSim = 100;
        int threshold = 25;

        //
        // mua3 = 1e-6
        foldername = "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Spacial_model/A3_mutation_rate_1e-6/Single_A3_Off/";

        sb = new StringBuilder();
        seedsFile = new FileIO(foldername + "_infos.csv", "w");

        for (int i = 0; i < nSim;i++){
            A3Parameters2D parameters = new A3Parameters2D();
            parameters.delete_thresh = threshold;

            int seed = SingleSim(parameters,totalTime,modifier,0,true,foldername);
            if (i==0){
                sb.append("br" + (int) (parameters.birth_rate*100) +
                        "dr" + (int) (parameters.death_rate*100) +
                        "a3l" + (int) (parameters.a3l*100 ) +
                        "a3m" + (int) (parameters.a3m*100 ) +
                        "a3h" + (int) (parameters.a3h*100 ) +
                        "seed" + seed +
                        "sim0" + "\n");
            }

            sb.append(seed + "\n");

        }
        seedsFile.Write(sb.toString());
        seedsFile.Close();

        //
        foldername = "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Spacial_model/A3_mutation_rate_1e-6/Single_A3_On/";

        sb = new StringBuilder();
        seedsFile = new FileIO(foldername + "_infos.csv", "w");

        for (int i = 0; i < nSim;i++){
            A3Parameters2D parameters = new A3Parameters2D();
            parameters.delete_thresh = threshold;
            parameters.delta = 1e-6;

            int seed = SingleSim(parameters,totalTime,modifier,0,true,foldername);
            if (i==0){
                sb.append("br" + (int) (parameters.birth_rate*100) +
                        "dr" + (int) (parameters.death_rate*100) +
                        "a3l" + (int) (parameters.a3l*100 ) +
                        "a3m" + (int) (parameters.a3m*100 ) +
                        "a3h" + (int) (parameters.a3h*100 ) +
                        "seed" + seed +
                        "sim0" + "\n");
            }

            sb.append(seed + "\n");

        }
        seedsFile.Write(sb.toString());
        seedsFile.Close();

        //
        foldername = "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Spacial_model/A3_mutation_rate_1e-6/Uniform_A3_On/";

        sb = new StringBuilder();
        seedsFile = new FileIO(foldername + "_infos.csv", "w");

        for (int i = 0; i < nSim;i++){
            A3Parameters2D parameters = new A3Parameters2D();
            parameters.delete_thresh = threshold;
            parameters.delta = 1e-6;
            parameters.a3l = 0.05;
            parameters.a3m = 0.9;
            parameters.a3h = 0.05;

            int seed = SingleSim(parameters,totalTime,modifier,0,true,foldername);
            if (i==0){
                sb.append("br" + (int) (parameters.birth_rate*100) +
                        "dr" + (int) (parameters.death_rate*100) +
                        "a3l" + (int) (parameters.a3l*100 ) +
                        "a3m" + (int) (parameters.a3m*100 ) +
                        "a3h" + (int) (parameters.a3h*100 ) +
                        "seed" + seed +
                        "sim0" + "\n");
            }

            sb.append(seed + "\n");

        }
        seedsFile.Write(sb.toString());
        seedsFile.Close();

        //
        foldername = "C:/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Spacial_model/A3_mutation_rate_1e-6/Skewed_A3_On/";

        sb = new StringBuilder();
        seedsFile = new FileIO(foldername + "_infos.csv", "w");

        for (int i = 0; i < nSim;i++){
            A3Parameters2D parameters = new A3Parameters2D();
            parameters.delete_thresh = threshold;
            parameters.delta = 1e-6;
            parameters.a3l = 0.3;
            parameters.a3m = 0.6;
            parameters.a3h = 0.1;

            int seed = SingleSim(parameters,totalTime,modifier,0,true,foldername);
            if (i==0){
                sb.append("br" + (int) (parameters.birth_rate*100) +
                        "dr" + (int) (parameters.death_rate*100) +
                        "a3l" + (int) (parameters.a3l*100 ) +
                        "a3m" + (int) (parameters.a3m*100 ) +
                        "a3h" + (int) (parameters.a3h*100 ) +
                        "seed" + seed +
                        "sim0" + "\n");
            }

            sb.append(seed + "\n");

        }
        seedsFile.Write(sb.toString());
        seedsFile.Close();

        return;
    }

    public static void DrawCells(A3MH2D model, UIGrid visCells) {
        // color half by drivers and half by passengers
        for (int i = 0; i < visCells.length; i++) {
            A3Cell2D c=model.GetAgent(i);

            if(c==null){
                visCells.SetPix(i, RGBA256(255,255,255, 255));
            } else{

                double knScaler = (model.KnMAX == 0) ? 255:255/(double)(model.KnMAX);
                double koScaler = (model.KoMAX == 0) ? 255:255/(double)(model.KoMAX);
                double kpScaler = (model.KpMAX == 0) ? 255:255/(double)(model.KpMAX);

                visCells.SetPix(i,()->{return RGBA256((int)(c.kp*kpScaler), (int)(c.ko*koScaler), (int)(c.kn*knScaler), 255);});
            }
        }
    }

    public static double GetNumGenotypes(A3MH2D model) {
        int theMAX = model.KpMAX*(model.KnMAX+1)*(model.KoMAX+1) + model.KnMAX*(model.KoMAX+1) + model.KoMAX;
        int kVec[] = new int[model.KpMAX*(model.KnMAX+1)*(model.KoMAX+1) + model.KnMAX*(model.KoMAX+1) + model.KoMAX+1];
        for (A3Cell2D c : model) {
            int ii = c.kp*(model.KnMAX+1)*(model.KoMAX+1) + c.kn*(model.KoMAX+1) + c.ko;
            kVec[ii]++;
        }

        double cnt = 0;
        for (int i = 0; i <= theMAX; i++) {
            if (kVec[i] > 0) {
                cnt++;
            }
        }
        return cnt;
    }

    public static double GetFitDiversity(A3MH2D model) {
        int theMAX = model.KpMAX*(model.KnMAX+1) + model.KnMAX;
        int kVec[] = new int[model.KpMAX*(model.KnMAX+1) + model.KnMAX + 1];
        for (A3Cell2D c : model) {
            int ii = c.kp*(model.KnMAX+1) + c.kn;
            kVec[ii]++;
        }

        double sum = 0;
        for (int i = 0; i <= theMAX; i++) {
            if (kVec[i] > 0) {
                sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
            }
        }
        return Math.exp(- sum);
    }

    public static double GetGenotypeDiversity(A3MH2D model) {
        int theMAX = model.KpMAX*(model.KnMAX+1)*(model.KoMAX+1) + model.KnMAX*(model.KoMAX+1) + model.KoMAX;
        int kVec[] = new int[model.KpMAX*(model.KnMAX+1)*(model.KoMAX+1) + model.KnMAX*(model.KoMAX+1) + model.KoMAX+1];
        for (A3Cell2D c : model) {
            int ii = c.kp*(model.KnMAX+1)*(model.KoMAX+1) + c.kn*(model.KoMAX+1) + c.ko;
            kVec[ii]++;
        }

        double sum = 0;
        for (int i = 0; i <= theMAX; i++) {
            if (kVec[i] > 0) {
                sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
            }
        }
        return Math.exp(- sum);
    }

    public static double GetDiversity(A3MH2D model, int passBool) {
        if (passBool == 0) {
            // count neg!
            int kVec[] = new int[model.KnMAX + 1];
            for (A3Cell2D c : model) {
                kVec[c.kn]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.KnMAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        } else if (passBool == 1){
            // count neu!
            int kVec[] = new int[model.KoMAX + 1];
            for (A3Cell2D c : model) {
                kVec[c.ko]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.KoMAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        } else if (passBool == 2){
            // count pos!
            int kVec[] = new int[model.KpMAX + 1];
            for (A3Cell2D c : model) {
                kVec[c.kp]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.KpMAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        } else if (passBool == 3){
            // count neg a3!
            int kVec[] = new int[model.Kna3MAX + 1];
            for (A3Cell2D c : model) {
                kVec[c.kna3]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.Kna3MAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        } else if (passBool == 4){
            // count neu a3!
            int kVec[] = new int[model.Koa3MAX + 1];
            for (A3Cell2D c : model) {
                kVec[c.koa3]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.Koa3MAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        } else if (passBool == 5){
            // count pos a3!
            int kVec[] = new int[model.Kpa3MAX + 1];
            for (A3Cell2D c : model) {
                kVec[c.kpa3]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.Kpa3MAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        } else if (passBool == 6){
            // count neg na3!
            int kVec[] = new int[model.Knna3MAX + 1];
            for (A3Cell2D c : model) {
                kVec[c.knna3]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.Knna3MAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);

        } else if (passBool == 7){
            // count neu na3!
            int kVec[] = new int[model.Kona3MAX + 1];
            for (A3Cell2D c : model) {
                kVec[c.kona3]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.Kona3MAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);
        } else {
            // count pos na3!
            int kVec[] = new int[model.Kpna3MAX + 1];
            for (A3Cell2D c : model) {
                kVec[c.kpna3]++;
            }

            // calculate diversity
            double sum = 0;
            for (int i = 0; i <= model.Kpna3MAX; i++) {
                if (kVec[i] > 0) {
                    sum += (double)((double)kVec[i] / (double)model.Pop())*Math.log10((double)kVec[i] / (double)model.Pop());
                }
            }
            return Math.exp(- sum);
        }
    }

    public static int CountMutators(A3MH2D model, int passBool) {
        int sum = 0;
        if (passBool == 0) {
            // count drivers!

            for (A3Cell2D c : model) {
                sum += c.kn;
            }
            return sum;

        } else if (passBool == 1) {
            // count passengers!
            for (A3Cell2D c : model) {
                sum += c.ko;
            }

            return sum;
        } else if (passBool == 2) {
            // count passengers!
            for (A3Cell2D c : model) {
                sum += c.kp;
            }

            return sum;
        } else if (passBool == 3) {
            // count passengers!
            for (A3Cell2D c : model) {
                sum += c.kna3;
            }

            return sum;
        } else if (passBool == 4) {
            // count passengers!
            for (A3Cell2D c : model) {
                sum += c.koa3;
            }

            return sum;
        } else if (passBool == 5) {
            // count passengers!
            for (A3Cell2D c : model) {
                sum += c.kpa3;
            }

            return sum;
        } else if (passBool == 6) {
            // count passengers!
            for (A3Cell2D c : model) {
                sum += c.knna3;
            }

            return sum;
        } else if (passBool == 7) {
            // count passengers!
            for (A3Cell2D c : model) {
                sum += c.kona3;
            }

            return sum;
        } else {
            // count passengers!
            for (A3Cell2D c : model) {
                sum += c.kpna3;
            }

            return sum;
        }
    }

    public static void BuildEverythingSingle(FileIO everythingDivsOutputFile, double[][] everythingDivs, int totalTime, int modifier) {
        StringBuilder sb = new StringBuilder();

        for (int jj = 0; jj < 32; jj++) {
            sb = new StringBuilder();
            for (int ii = 0; ii < (totalTime / modifier); ii++) {
                sb.append(everythingDivs[jj][ii] + ",");
            }
            sb.append(everythingDivs[jj][(totalTime / modifier)] + "\n");
            everythingDivsOutputFile.Write(sb.toString());
        }
        everythingDivsOutputFile.Close();
    }

    public static void SaveAll(FileIO parentsFile, FileIO geneticMullerFile, FileIO allStatusFile,
                               FileIO negStatusFile, FileIO neuStatusFile, FileIO posStatusFile,
                               FileIO nega3StatusFile, FileIO neua3StatusFile, FileIO posa3StatusFile,
                               FileIO negna3StatusFile, FileIO neuna3StatusFile, FileIO posna3StatusFile,
                               FileIO genotypeStatusFile, FileIO fitnessStatusFile,
                               int[][] mullerGenetic, int[] progenyToParentIDs, int[] all_status,
                               int[] neg_status, int[] neu_status, int[] pos_status,
                               int[] nega3_status, int[] neua3_status, int[] posa3_status,
                               int[] negna3_status, int[] neuna3_status, int[] posna3_status,
                               int[] genotype_status, int[] fitness_status,
                               int expectedProgeny, int maxProgenyID, int totalTime, int modifier,
                               int threshold) {

        StringBuilder sb;
        int[] maxNumClones = new int[maxProgenyID];
        boolean[] cloneIDtoSave = new boolean[maxProgenyID];
        cloneIDtoSave[0] = false;
        cloneIDtoSave[1] = true;

        int maxNum;
        for (int jj = 1; jj < maxProgenyID; jj++){
            maxNum = 0;
            for (int ii = 0; ii < (totalTime / modifier)+1; ii++) {
                maxNum = (mullerGenetic[jj][ii] > maxNum) ? mullerGenetic[jj][ii]: maxNum;
            }
            maxNumClones[jj] = maxNum;

            if (maxNum >= threshold) { cloneIDtoSave[jj] = true; }

        }

        for (int jj = maxProgenyID-1; jj > 1; jj--){
            if (cloneIDtoSave[jj]){
                cloneIDtoSave[progenyToParentIDs[jj]] = true;
            }
        }


        //
        sb = new StringBuilder();

        for (int time = 0; time < (totalTime / modifier); time++) { sb.append(time*modifier + ","); }
        sb.append(modifier*((totalTime / modifier)) + "\n");
        geneticMullerFile.Write(sb.toString());
        for (int jj = 1; jj < maxProgenyID; jj++) {
            if (cloneIDtoSave[jj]) {
                sb = new StringBuilder();
                sb.append(jj + ",");
                for (int ii = 0; ii < (totalTime / modifier); ii++) {
                    sb.append(mullerGenetic[jj][ii] + ",");
                }
                sb.append(mullerGenetic[jj][(totalTime / modifier)] + "\n");
                geneticMullerFile.Write(sb.toString());
            }
        }
        geneticMullerFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(progenyToParentIDs[jj] + "\n");
            }
        }
        parentsFile.Write(sb.toString());
        parentsFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(all_status[jj] + "\n");
            }
        }
        allStatusFile.Write(sb.toString());
        allStatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(neg_status[jj] + "\n");
            }

        }
        negStatusFile.Write(sb.toString());
        negStatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(neu_status[jj] + "\n");
            }

        }
        neuStatusFile.Write(sb.toString());
        neuStatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(pos_status[jj] + "\n");
            }

        }
        posStatusFile.Write(sb.toString());
        posStatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(nega3_status[jj] + "\n");
            }

        }
        nega3StatusFile.Write(sb.toString());
        nega3StatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(neua3_status[jj] + "\n");
            }

        }
        neua3StatusFile.Write(sb.toString());
        neua3StatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(posa3_status[jj] + "\n");
            }

        }
        posa3StatusFile.Write(sb.toString());
        posa3StatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(negna3_status[jj] + "\n");
            }

        }
        negna3StatusFile.Write(sb.toString());
        negna3StatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(neuna3_status[jj] + "\n");
            }

        }
        neuna3StatusFile.Write(sb.toString());
        neuna3StatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(posna3_status[jj] + "\n");
            }

        }
        posna3StatusFile.Write(sb.toString());
        posna3StatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(genotype_status[jj] + "\n");
            }

        }
        genotypeStatusFile.Write(sb.toString());
        genotypeStatusFile.Close();

        //
        sb = new StringBuilder();
        for (int jj = 1; jj < maxProgenyID; jj++){
            if (cloneIDtoSave[jj]) {
                sb.append(fitness_status[jj] + "\n");
            }

        }
        fitnessStatusFile.Write(sb.toString());
        fitnessStatusFile.Close();

    }

}