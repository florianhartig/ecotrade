unit Unit1mitpreis;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  ExtCtrls, clipbrd, jpeg, StdCtrls, filectrl, inifiles, ComCtrls, math,
  TeEngine, Series, TeeProcs, Chart, Menus;

type
  TForm1 = class(TForm)
    Image1: TImage;
    GroupBoxEcologicalModel: TGroupBox;
    PageControlStatistics: TPageControl;
    TabSheet1: TTabSheet;
    TabSheet2: TTabSheet;
    TabSheet3: TTabSheet;
    GroupBoxRunSettings: TGroupBox;
    EditPrecision: TEdit;
    Label28: TLabel;
    Edit10: TEdit;
    Label10: TLabel;
    Edit3: TEdit;
    Label3: TLabel;
    CheckBoxShowWarnings: TCheckBox;
    GroupBoxVariation: TGroupBox;
    Label19: TLabel;
    Label16: TLabel;
    Label17: TLabel;
    Label18: TLabel;
    Label24: TLabel;
    varyRun: TButton;
    Filename: TEdit;
    vAround: TEdit;
    vRange: TEdit;
    vSteps: TEdit;
    variationTries: TEdit;
    varound2d1: TEdit;
    vRange2d1: TEdit;
    vSteps2d1: TEdit;
    varound2d2: TEdit;
    vRange2d2: TEdit;
    vsteps2d2: TEdit;
    RunStatistics: TListBox;
    TabSheet4: TTabSheet;
    TabSheet5: TTabSheet;
    Neighboursshow: TListBox;
    RadioGroupStartDistribution: TRadioGroup;
    GroupBoxControl: TGroupBox;
    Button2: TButton;
    EditInifilename: TEdit;
    Button3: TButton;
    ButtonSaveIni: TButton;
    GroupBoxEconomics: TGroupBox;
    Edit1: TEdit;
    Label1: TLabel;
    RadioGroupCosts: TRadioGroup;
    Label2: TLabel;
    Edit2: TEdit;
    Label11: TLabel;
    Edit11: TEdit;
    Label7: TLabel;
    Edit7: TEdit;
    GroupBoxStatistics: TGroupBox;
    Label4: TLabel;
    Edit4: TEdit;
    Label14: TLabel;
    goaldiff: TEdit;
    CostLevel: TEdit;
    Label23: TLabel;
    Edit8: TEdit;
    Label8: TLabel;
    Costs: TEdit;
    Label13: TLabel;
    EditAveragePayments: TEdit;
    Label21: TLabel;
    Edit5: TEdit;
    Label5: TLabel;
    Edit6: TEdit;
    Label6: TLabel;
    Turnover: TEdit;
    Label12: TLabel;
    EditMortality: TEdit;
    EditMortalityCorrelation: TEdit;
    EditDispersalRate: TEdit;
    EditDispersalDistance: TEdit;
    Label15: TLabel;
    Label22: TLabel;
    Label27: TLabel;
    Label29: TLabel;
    ButtonOpenIni: TButton;
    GroupBoxGraphicWindowControlls: TGroupBox;
    EditZoom: TEdit;
    Label30: TLabel;
    ShowPopulatedHabitats: TEdit;
    Button1: TButton;
    CheckBoxEconomicModelOn: TCheckBox;
    EditIntensity: TEdit;
    CheckBoxEcologicalModelOn: TCheckBox;
    ButtonEcologicalModelReset: TButton;
    ShowColonizationProbability: TEdit;
    Label31: TLabel;
    Label32: TLabel;
    CheckBoxShowGrafics: TCheckBox;
    UpDownZoom: TUpDown;
    RadioGroupSerial: TRadioGroup;
    ChartExtinction: TChart;
    RadioGroupShowValues: TRadioGroup;
    Series1: TFastLineSeries;
    EditAltruism: TEdit;
    CheckBoxAltruismon: TCheckBox;
    Label33: TLabel;
    Series2: TFastLineSeries;
    MainMenu1: TMainMenu;
    test1: TMenuItem;
    rt1: TMenuItem;
    Test21: TMenuItem;
    PopupMenuVariation1: TPopupMenu;
    AgglomerationBonusw1: TMenuItem;
    CostVariation1: TMenuItem;
    Dimension1: TMenuItem;
    Stepsbetweenupdating1: TMenuItem;
    GroupBoxOptimization: TGroupBox;
    CheckBoxOptimizationOn: TCheckBox;
    EditOptimizationPercentage: TEdit;
    RadioGroupTradeConstraints: TRadioGroup;
    CheckBoxIndependentrunon: TCheckBox;
    N1: TMenuItem;
    PopupMenuVariation2: TPopupMenu;
    Popupmenuvariation2item0: TMenuItem;
    TabSheet6: TTabSheet;
    ButtonGTRun: TButton;
    ButtonGTReset: TButton;
    EditGTcc: TEdit;
    EditGTdc: TEdit;
    EditGTcd: TEdit;
    EditGTdd: TEdit;
    ShowSDcosts: TEdit;
    Showmcostchange: TEdit;
    Label25: TLabel;
    Label26: TLabel;
    EditCostSize: TEdit;
    Label34: TLabel;
    Editmradius: TEdit;
    Editvradius: TEdit;
    Label35: TLabel;
    EditDispersalPreference: TEdit;
    Label36: TLabel;
    Chart1: TChart;
    Series3: TFastLineSeries;
    GroupBoxEconomicCharts: TGroupBox;
    EditCostDistributionRange: TEdit;
    CheckBoxCalculateCostDistribution: TCheckBox;
    Series4: TFastLineSeries;
    ButtonIniSaveAs: TButton;
    Label37: TLabel;
    CheckBoxSaveVideo: TCheckBox;
    EditStartTemperature: TEdit;
    EditTemperatureDecay: TEdit;
    Label38: TLabel;
    Label39: TLabel;
    CheckBoxOptimizationShow: TCheckBox;
    Label20: TLabel;
    TradeRuns: TEdit;
    RadioGroupVariationChoice: TRadioGroup;
    ButtonBreak: TButton;
    ButtonLogWindowClear: TButton;
    TabSheetOptimalSpatialBonus: TTabSheet;
    Chart2: TChart;
    Series5: TBarSeries;
    LabelVariationParameter: TLabel;
    Label40: TLabel;
    EditLength: TEdit;
    Edit9: TEdit;
    Label9: TLabel;
    LabelNeighbors: TLabel;
    Memo1: TMemo;
    GroupBoxEcologicalChartControls: TGroupBox;
    CheckBoxShowMortality: TCheckBox;
    CheckBoxShowPopulation: TCheckBox;
    CheckBoxShowNewlyCreatedHabitats: TCheckBox;
    CheckBoxShowDeaths: TCheckBox;
    CheckBoxShowDestroyedHabitats: TCheckBox;
    Button4: TButton;
    Button5: TButton;
    Label41: TLabel;
    EditReduce: TEdit;
    LabelRed: TLabel;
    Button6: TButton;
    procedure Edit1Change(Sender: TObject);
    procedure Edit2Change(Sender: TObject);
    procedure Edit3Change(Sender: TObject);
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure Edit7Change(Sender: TObject);
    procedure Edit9Change(Sender: TObject);
    procedure Edit10Change(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure Edit11Change(Sender: TObject);
    procedure TradeRunsChange(Sender: TObject);
    procedure varyRunClick(Sender: TObject);
    procedure vRangeChange(Sender: TObject);
    procedure vStepsChange(Sender: TObject);
    procedure variationTriesChange(Sender: TObject);
    procedure vRange2d1Change(Sender: TObject);
    procedure vSteps2d1Change(Sender: TObject);
    procedure vRange2d2Change(Sender: TObject);
    procedure vSteps2d2Change(Sender: TObject);
    procedure FilenameChange(Sender: TObject);
    procedure RadioGroupShowValuesClick(Sender: TObject);
    procedure EditIntensityChange(Sender: TObject);
    procedure EditPrecisionChange(Sender: TObject);
    procedure ButtonSaveIniClick(Sender: TObject);
    procedure CheckBoxShowWarningsClick(Sender: TObject);
    procedure form1close(Sender: TObject; var Action: TCloseAction);
    procedure EditMortalityChange(Sender: TObject);
    procedure EditMortalityCorrelationChange(Sender: TObject);
    procedure EditDispersalDistanceChange(Sender: TObject);
    procedure EditDispersalRateChange(Sender: TObject);
    procedure ButtonOpenIniClick(Sender: TObject);
    procedure CheckBoxEcologicalModelOnClick(Sender: TObject);
    procedure EditZoomChange(Sender: TObject);
    procedure ButtonEcologicalModelResetClick(Sender: TObject);
    procedure CheckBoxEconomicModelOnClick(Sender: TObject);
    procedure UpDownZoomClick(Sender: TObject; Button: TUDBtnType);
    procedure RadioGroupSerialClick(Sender: TObject);
    procedure RadioGroupCostsClick(Sender: TObject);
    procedure RadioGroupStartDistributionClick(Sender: TObject);
    procedure EditAltruismChange(Sender: TObject);
    procedure CheckBoxAltruismonClick(Sender: TObject);
    procedure CheckBoxOptimizationOnClick(Sender: TObject);
    procedure EditOptimizationPercentageChange(Sender: TObject);
    procedure RadioGrouptradeconstraintsClick(Sender: TObject);
    procedure CheckBoxIndependentrunonClick(Sender: TObject);
    procedure N1stvariationClick(Sender: TObject);
    procedure vAroundChange(Sender: TObject);
    procedure Popupmenuvariation2itemClick(Sender: TObject);
    procedure ButtonGTRunClick(Sender: TObject);
    procedure ButtonGTResetClick(Sender: TObject);
    procedure EditGTccChange(Sender: TObject);
    procedure EditGTcdChange(Sender: TObject);
    procedure EditGTdcChange(Sender: TObject);
    procedure EditGTddChange(Sender: TObject);
    procedure EditCostSizeChange(Sender: TObject);
    procedure EditmradiusChange(Sender: TObject);
    procedure EditvradiusChange(Sender: TObject);
    procedure EditDispersalPreferenceChange(Sender: TObject);
    procedure Chart1Click(Sender: TObject);
    procedure Image1Click(Sender: TObject);
    procedure ButtonIniSaveAsClick(Sender: TObject);
    procedure EditStartTemperatureChange(Sender: TObject);
    procedure EditTemperatureDecayChange(Sender: TObject);
    procedure CheckBoxOptimizationShowClick(Sender: TObject);
    procedure varound2d1Change(Sender: TObject);
    procedure varound2d2Change(Sender: TObject);
    procedure RadioGroupVariationChoiceClick(Sender: TObject);
    procedure ButtonBreakClick(Sender: TObject);
    procedure ButtonLogWindowClearClick(Sender: TObject);
    procedure EditLengthChange(Sender: TObject);
    procedure Memo1Change(Sender: TObject);
    procedure Button4Click(Sender: TObject);
    procedure Button5Click(Sender: TObject);
    procedure EditReduceChange(Sender: TObject);
    procedure Button6Click(Sender: TObject);


  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;

implementation

uses Unit2 ;

{$R *.DFM}



const maxdim  = 50;
type realmatrix=array[0..maxdim+1,0..maxdim+1] of double;
     intmatrix=array[0..maxdim+1,0..maxdim+1] of integer;
     dynrealmatrix=array of array of double;
     coordinate = array[0..1] of integer;
     neighborcoordinates = array [0..8] of coordinate;
     coordinatematrix = array [0..maxdim*maxdim] of coordinate;
     coordinatepointer = ^coordinate;
     Pintmatrix = ^intmatrix;


// VARIABLES --------------------------------------------------------------------------

var //state variables ----------------------------------------

    c,cneighbour: realmatrix;       // costs for each patch
    m, mtemp: intmatrix;            // number of neighbours for each patch
    x, xtemp, xold : intmatrix;           // habitats matrices
    p : intmatrix;
    ptemp : realmatrix;
    habitats, populatedhabitats, unpopulatedhabitats, nonhabitats : coordinatematrix;

    fitnessmatrix : realmatrix;           // matrix for the game theory tryout
    mortalitymatrix : realmatrix;         // for ecomodel
    deathmatrix : intmatrix;              // for ecomodel
    strategymatrix : realmatrix;          // for testin the best evolutionary strategie

                                          // holds the neighbors
    neighbormatrix : array of array of integer ;

    //statistics -----------------------------------------

    price: double ;         //price per patch
    xtot : integer;                 //total amount of patches
    totalcosts : double ;
    turnoverrate : double;
    calculations : boolean ;
    mprecision: double;


    // m is for mean and is an aggregated value
    mxtot, mprice, meantotalcosts : double;  // aggregated values in one run


    averagecostlevel : double;
    averageneigbours : double ;
    sdcosts,mcostchange : double; // calculates the sd and the mean change of the costs


    neighbourstatistics : array of double; // percentage of neighbours
    ztot: double;
    loops, loopsum : integer; // number of loops run to get the price, just for performance controll
    outfile: TextFile;
    count,turncount:integer;
    summary  : array[1..10,1..2] of double;
    connectivity : double;
    numberofhabitats : integer ;
    numberofpopulatedhabitats: integer = 0;numberofunpopulatedhabitats : integer ;
    disppref : array of double; // for storing the disperal preferences in dispers
    mcolonizationsum : double;
    storetime: Tdatetime;
    mortalitystatistics, populationstatistics : array of double;
    ecocounter : integer;
    neighbornumber : integer;


    // variation ----------------

    variation1,variation2 : integer;

    varyaround : double;
    varyRange : double;            //range for variation
    varysteps : integer ;          // variation steps
    varytries : integer ;          //number of independent tries in the variation

    varyaround2d1 : double;
    varyRange2d1 : double;         //range for variation
    varysteps2d1 : integer ;       // variation steps
    varyaround2d2 : double;
    varyRange2d2 : double;         //range for variation
    varysteps2d2 : integer ;       // variation steps

    variationChoice : integer ;
    variationParameterChoice : integer ;
    reduce : integer;              // for running several population counts after the burnin of the economic model


    // global and GUI options

    showvalues : boolean;
    showwarnings : boolean ;
    ecologicalmodelon: boolean;
    economicmodelon: boolean;
    statisticsfilename : string;
    inifilename : string ;
    remarks : string;
    altruismon : boolean ;
    independentrunon : boolean;
    circularx : boolean;
    circulary : boolean;
    breakrun : boolean;

    // graphics ------------------

    bild: TBitmap;                 // for the graphical output
    zoom: integer;                 // size of the picture


    // Parameters for the economic model
    //----------------------------

    lambda: double ;               //desired level of habitats
    dim: integer ;                 //dimension of the patch matrix
    w: double ;                    //interaction strength
    interactionlength: double ;    // length in the connectivity measures
    trades : integer;              // number of trades without variation of costs
    maxt: integer ;                // number of runs inside the statistics calculation
    runs: integer ;                // number of runs outside the statistics calculations
    sigma: double ;                //distribution scale of costs
    pricestability: double ;       //controlls the range of the prices
    intensity : integer;           //fraction of patches traded per round


    // the following are calculated on runtime on base of the parameters
    desiredBenefit: double;  sqrdim:integer;
    benefitarray : array of double;
    distancematrix : array[0..maxdim,0..maxdim] of double;


    // Switches of the economic model

    serial: integer;               // serial or parallel updating in trade
    precision: double;             // precision for the meeting of the constraints
    costvariation:integer;         // switch for the generation of cost surfaces
    startdistribution :integer;    // initialization
    altruism : double;             // coordination of neighbors

    tradeconstraints : integer;    // switch between fixing budget or target
    costsize : integer;

    //-----------------------------------------------
    // ecological model

    mortality : double ;           //mortality
    mradius,vradius: double ;      // for circles
    rho : double;                  //mortality correlation
    dispersalrate : double ;       //dispersal rate
    alpha : double;                //dispersal distance
    dispersalradius : double;      // radius for dispersal
    dispersalpreference : double;


    // --------------------------------------------------
    // game theory

     cc,cd,dc,dd : double;

     //---------------------------------------------------
     // optimization

     optimizationon:boolean;
     optimizationpercentage:double ;
     optimizationshow : boolean;
     optimizationAlgorithm : integer;
     StartTemperature, TemperatureDecay : double;      //for SA algorithm


// forward declarations

procedure enableforms; forward;
procedure warning(code:integer); forward;
function dispersalProbability (inp:double):double; forward;
procedure write_ini(filename:string); forward;
procedure read_ini  ; forward;
procedure reset ; forward;
procedure run(cos:integer); forward;




//----------------------------------------------------------
//help procedures-----------------------------------------

function sign (inp:double):integer ;
begin
if inp > 0 then result := +1
else result := -1;
end;


function modx(i:integer):integer;
begin
  result := (i+dim-1) mod dim + 1 ;
end;

function neighbours(i,j : integer): neighborcoordinates;
// calcuates the neighbors of a coordinate, result is an array with neighborcoordinates
var left,right,up,down : integer;
begin
  left := (i+dim-2) mod dim +1 ;
  right :=  i mod dim +1 ;
  up  :=  (j+dim-2) mod dim +1 ;
  down  := j mod dim +1;
  result[0,0]:=i             ;result[0,1]:= j  ;
  result[1,0]:=right         ;result[1,1]:= j;
  result[2,0]:=right         ;result[2,1]:= down;
  result[3,0]:=i             ;result[3,1]:= down ;
  result[4,0]:=left          ;result[4,1]:= down ;
  result[5,0]:=left          ;result[5,1]:=  j ;
  result[6,0]:=left          ;result[6,1]:=  up;
  result[7,0]:=i             ;result[7,1]:=  up;
  result[8,0]:=right         ;result[8,1]:=  up ;
end;

function sqrdistance(i,j,k,l : integer):double;  // euklidian distance on the grid, periodic boundary conditions
var xx,yy : integer;
begin
  if abs (i-k) > dim/2 then xx := abs (abs(i-k)-dim) else xx := abs (i-k);
  if abs (j-l) > dim/2 then yy := abs (abs(j-l)-dim) else yy := abs (j-l);
  result := sqr (xx) + sqr (yy);
end;

function distance(i,j,k,l : integer):double;    // euklidian distance on the grid, periodic boundary conditions
var xx,yy : integer;
begin
  if abs (i-k) > dim/2 then xx := abs (abs(i-k)-dim) else xx := abs (i-k);
  if abs (j-l) > dim/2 then yy := abs (abs(j-l)-dim) else yy := abs (j-l);
  result := sqrt( sqr (xx) + sqr (yy));
end;

function idToCoordinate(i: integer):coordinate;  // coordinate id counts from 1 to sqrdim
begin
  result[0] := (i) mod dim +1 ;
  result[1] := (i) div dim +1 ;
end;

function coordinatToId(i,j:integer):integer;
begin
  result := (j-1)*dim + (i -1 )
end;

procedure setneighbors;   // calculates the matrix with the neighbors of each cell
var i,j, count: integer;
    tempc, tempn : coordinate;
        dist : double;
begin
   neighbornumber:= 0;
   for i:= 1 to sqrdim-1 do
   begin
    tempc := idToCoordinate(i);
    dist := distance(1, 1, tempc[0], tempc[1]);
    if  dist <= interactionlength then neighbornumber := neighbornumber + 1;
   end;
   SetLength(neighbormatrix, sqrdim, neighbornumber);
   for i:= 0 to sqrdim-1 do
   begin
     tempc := idToCoordinate(i);
     count := 0;
     for j:= 0 to sqrdim-1 do
     begin
       if (i <> j) then
       begin
         tempn := idToCoordinate(j);
         dist := distance(tempn[0], tempn[1], tempc[0], tempc[1]);
         if  dist <= interactionlength then
         begin
           neighbormatrix[i,count] := j;
           count := count +1 ;
         end;
       end;
     end;
  end;
  setlength(neighbourstatistics,neighbornumber+1);
end;

procedure calculateNumberOfNeighbors_x_m  ;
//calculates the interaction benefit for every patch on x
var i,k, counter: integer;
     tempc : coordinate;
begin
   for i := 0 to sqrdim-1 do
   begin
     counter := 0;
     for k:= 0 to neighbornumber-1 do
     begin
      tempc := idToCoordinate(neighbormatrix[i,k]);
      counter := counter + x[tempc[0], tempc[1]];
     end;
     tempc := idToCoordinate(i);
     m[tempc[0], tempc[1]] := counter;
   end;
end;

procedure calculateNumberOfNeighbors_xtemp_mtemp ;
//calculates the interaction benefit for every patch  on xtemp
var i,k, counter: integer;
     tempc : coordinate;
begin
   for i := 0 to sqrdim-1 do
   begin
     counter := 0;
     for k:= 0 to neighbornumber-1 do
     begin
      tempc := idToCoordinate(neighbormatrix[i,k]);
      counter := counter + xtemp[tempc[0], tempc[1]];
     end;
     tempc := idToCoordinate(i);
     mtemp[tempc[0], tempc[1]] := counter;
   end;
end;


procedure calculatesinglemtemp(k,l:integer);
// calculates the changes of the mtemp matrix after a change of a singe xtemp of which the
// coordinates are passed to the procedure
var m,n,i,j : integer;
var co, co2 : neighborcoordinates;
begin
  co := neighbours(k,l);
  for m:= 1 to 8 do
  begin
    i := co[m,0]; j:= co[m,1];
    co2 := neighbours(i,j);
    mtemp[i,j]:=0;
    for n:= 1 to 8 do mtemp[i,j]:=mtemp[i,j]+xtemp[co2[n,0],co2[n,1]];
  end;
end;


// --------------------------------------------------------------
// grid procedures

procedure randomize_costs(par:integer) ;
//is called with an integer which decides how the cost variation is done
// 0 = assigns randomn costs between 1-sigma and 1+sigma to the patches;
// 1 = random walk - pricestability*(1-c[i,j]) ;
// 2 = random walk with interaction
// 3 = no variation
// 4
// 5  ... see below
// 6

var i,j,ii,jj : integer;
    part: double;
    lognormalmu,lognormalsigma : double;
begin
  case par of
    0:  //assigns randomn costs between 1-sigma and 1+sigma to the patches;
    begin
      for i:=1 to dim do for j:=1 to dim do c[i,j]:=1+sigma*(1-2*random);
    end;
    1:     //randomwalk
    begin
      for i:=1 to dim do for j:=1 to dim do
      begin
        //c[i,j]:= c[i,j]*(1/averagecostlevel)+ sigma*(1-2*random) + pricestability*(1-c[i,j]) ;
        c[i,j]:= c[i,j] + pricestability *(1-2*random)+ ( sigma  * sign(1 - c[i,j]) * Power( abs(1 - c[i,j]),0.5)) ;

        //c[i,j]:= c[i,j]*(1/averagecostlevel)+ sigma*(1-2*random) + 0.1*sigma*(1-c[i,j]) ;

        // if (  c[i,j]) > (1+sigma) then c[i,j] := 1+ sigma else
        // if ( c[i,j]) < (1-sigma) then c[i,j] := 1 - sigma;
      end;
    end;
    2:     // random walk with interaction
    begin
      for i:=1 to dim do for j:=1 to dim do cneighbour[i,j]:=0;
      for i:=1 to dim do for j:=1 to dim do
      begin
        if i-1>=1 then ii:=i-1 else ii:=dim; cneighbour[i,j]:=cneighbour[i,j]+c[ii,j]; cneighbour[ii,j]:=cneighbour[ii,j]+c[i,j];
        if j-1>=1 then jj:=j-1 else jj:=dim; cneighbour[i,j]:=cneighbour[i,j]+c[ii,jj]; cneighbour[ii,jj]:=cneighbour[ii,jj]+c[i,j];
        if j+1<=dim then jj:=j+1 else jj:=1; cneighbour[i,j]:=cneighbour[i,j]+c[ii,jj]; cneighbour[ii,jj]:=cneighbour[ii,jj]+c[i,j];
        if j-1>=1 then jj:=j-1 else jj:=dim; cneighbour[i,j]:=cneighbour[i,j]+c[i,jj];  cneighbour[i,jj]:=cneighbour[i,jj]+c[i,j];
      end ;
      for i:=1 to dim do for j:=1 to dim do
      begin
      c[i,j]:= ((1-pricestability)*c[i,j] + pricestability*cneighbour[i,j]/8)*(1/averagecostlevel) + randg(0,sigma) ;
      end;
    end;
    3: ;
    4:
    begin
      i:=1;j:=1;
      while i < (dim +1 ) do
      begin
        while j < (dim +1 )do
        begin
          part := 1+sigma*(1-2*random);
          for ii := i to i+ (costsize-1) do for jj := j to j + (costsize-1) do
          begin
            c[ii,jj]:= part;
          end;
          j := j + costsize;
        end;
      i := i + costsize;
      j:= 1;
      end;
    end;
    5:  //lognormal  distribution;
    begin
      lognormalmu:= ln(1)- 0.5* ln(1+sigma/1);
      lognormalsigma := sqrt(ln(1+sigma/1));
      begin
        for i:=1 to dim do for j:=1 to dim do c[i,j]:= exp(randg(lognormalmu, lognormalsigma));
      end;
    end ;
    6:  //gauss  distribution;
    begin
      begin
        for i:=1 to dim do for j:=1 to dim do c[i,j]:= randg(1, sigma);
      end;
    end ;
    else if showwarnings then form1.RunStatistics.Items.insert(0,'falscher Aufruf von randomize costs');
  end;
end;

procedure createHabitatList;
// create a lists with the coordinates of all conserved cells
var i,j: integer;
begin
  numberofhabitats:=0;
  for i:=1 to dim do  for j:=1 to dim do
  begin
    if x[i,j]=1 then
    begin
      habitats[numberofhabitats,0]:= i;
      habitats[numberofhabitats,1]:= j;
      numberofhabitats:= numberofhabitats + 1;
    end;
  end;
  if showwarnings and ( numberofhabitats <> xtot) then
  begin
      form1.RunStatistics.Items.insert(0,'Error in createHabitatList, number of habitats <> xtot');
  end;
end;

procedure createnonHabitatList;
// create a lists with the coordinates of all nonconserved cells
var i,j,k: integer;
begin
  k:=0;
  for i:=1 to dim do  for j:=1 to dim do
  begin
    if x[i,j]=0 then
    begin
      nonhabitats[k,0]:= i;
      nonhabitats[k,1]:= j;
      k:= k + 1;
    end;
  end;
end;


function ecological_benefit (inp:integer):double;
// calculates the ecological benefit ... it seems the function is way faster than
// an array with stored return values, at least for one floating point multiplication
begin
  // if inp > 5 then result := 2*w  else if inp > 3 then result := 8*w else result := 0;
  //if inp > 5 then result := 8*w else result := 0;
  result := (1-w) + w/neighbornumber * inp;
  //result := w*sqrt(inp)
end;

function ecological_benefit_alt (inp:integer):double; // altruistic version
// calculates the ecological benefit ... it seems the function is way faster than
// an array with stored return values, at least for one floating point multiplication
begin
        result := (1-w) + (1+altruism)*w/neighbornumber * inp;
end;

function marginal_ecological_benefit(inp:integer):double; // marginal version
// calculates the ecological benefit ... it seems the function is way faster than
// an array with stored return values, at least for one floating point multiplication
begin
        result := (1-w) + 2*w/neighbornumber * inp;
end;

{function ecological_benefite (inp:integer):double;
begin
  case inp of
  0..8 : result := benefitarray[inp]  ;
  else if showwarnings then showmessage('Error in ecological benefit') ;
  end;
end; }


procedure constants_calculate;
// calculates some standard values used later
var i,j : integer;
begin
  sqrdim := sqr(dim);
  setneighbors; // calculates the matrix of neighbors

  {for i := 0 to 8 do
  begin
    benefitarray[i]:= ecological_benefit(i)
  end; }
  desiredBenefit := lambda*(ecological_benefit(neighbornumber));
  //setlength(distancematrix, halfdim, halfdim);
  for i := 0 to dim-1 do for j := 0 to dim-1 do
  begin
    distancematrix[i,j]:= dispersalProbability(sqrt(i*i+j*j))
  end;
  if maxt < 1 then maxt := 1 ; // avoid floating point error for inside 
end;

procedure basic_calculations_do;
// basic calculations, does include calculations used at runtime
// but no calculations for the output
var i,j: integer;
begin
  calculateNumberOfNeighbors_x_m;
  xtot:=0; ztot:=0;averagecostlevel:=0;
  for i:=1 to dim do for j:=1 to dim do
  begin
    ztot:= ztot +  x[i,j]*(ecological_benefit(m[i,j]));
    xtot:=xtot+x[i,j];
    averagecostlevel := averagecostlevel + c[i,j] ;
  end;
  if xtot <= 0 then xtot :=1;
  ztot := ztot/sqrdim;
  averagecostlevel := averagecostlevel/sqrdim;
end;

procedure calculations_do;
//  calculates runtime and output values
var i,j: integer;
    neighbstat : array of integer;
    connect : integer;
    help : double;


begin
  // initializations ------------------------
  calculateNumberOfNeighbors_x_m;
  totalcosts:=0; xtot:=0; ztot:=0;
  setlength(neighbstat, neighbornumber+1);
  for i:=0 to neighbornumber do neighbstat [i]:=0;
  averagecostlevel:=0;
  connect := 0;
  averageneigbours := 0;
  sdcosts:= 0;
  mcostchange:= 0;
  help := 1  ;
  // calculations ---------------------
  for i:=1 to dim do for j:=1 to dim do
  begin
    ztot:= ztot +  x[i,j]*(ecological_benefit(m[i,j]));
    xtot:=xtot+x[i,j];
    totalcosts := totalcosts + x[i,j]*c[i,j];
    averagecostlevel := averagecostlevel + c[i,j] ;
    sdcosts:= sdcosts + abs(c[i,j]-1);
    mcostchange := mcostchange + abs(c[i,j]-help); help :=c[i,j];
    if x[i,j] = 1 then
    begin
    connect := connect + m[i,j] ;
    neighbstat[m[i,j]]:=neighbstat[m[i,j]]+1;
    end;
  end;
  if xtot <= 0 then xtot := 1;
  ztot := ztot/sqrdim;
  connectivity := connect/xtot/neighbornumber ;
  totalcosts := totalcosts/sqrdim;
  averagecostlevel := averagecostlevel/sqrdim;
  sdcosts := sdcosts/sqrdim;
  mcostchange := mcostchange/sqrdim;
  for i:=0 to neighbornumber do
  begin
    averageneigbours := averageneigbours + neighbstat[i] * (i-1) /xtot;
    neighbourstatistics [i]:=neighbstat[i]/xtot*100;
  end;
end;

procedure randomize_patches(par:integer) ;
//is called with an integer which decides how the initial patch distribution is done
// 0 = creates a randomn distribution of 0 and 1 in the patches
// 1 = creates random distribution of patches and half half distribution of costs
// 2 = creates half half distribution of the patches
// ... for more options see below
var i,j,ii,jj,iii,jjj,k: integer;
    radius, sqrradius : double;
begin
  case par of
    0:     // random
      for i:=1 to dim do for j:=1 to dim do if random < lambda then x[i,j]:= 1 else x[i,j]:= 0 ;
    1:     // parted costs
      begin
        for i:=1 to dim do for j:=1 to dim do
        begin
          x[i,j]:=round(random) ;
          if j>round(dim/2) then c[i,j]:= 1.5 else c[i,j]:= 0.5;
        end
      end;
    2:   // parted land
      begin
       for i:=1 to dim do for j:=1 to dim do
       if j>round(dim/2) then x[i,j]:= 1  else x[i,j]:= 0;
      end;
    3:    // all habitat
      for i:=1 to dim do for j:=1 to dim do
        x[i,j]:=1 ;
    4:      // gradient
      begin
        for i := 1 to dim  do for j:=1 to dim do c[i,j] := 1 - 2*abs(dim/2 - j)/dim;
        for i:=1 to dim do for j:=1 to dim do if random < lambda then x[i,j]:= 1 else x[i,j]:= 0 ;
       // constants_calculate;
       // basic_calculations_do;   count :=0;
        //run(costvariation);
       // calculations_do;
      end;
    5:    ;  // inifile
    6:      // circle
      begin
        for i:=1 to dim do for j:=1 to dim do x[i,j]:= 0;
        i := round(dim/2) ; j := round (dim/2);
        radius := sqrt(sqrdim* lambda /pi);
        sqrradius := sqr(radius);
        k:=0;
        for ii:=i-round(radius) to i + round(radius) do for jj:=j-round(radius) to j + round(radius) do
        begin
           if (sqr(ii-i) + sqr(jj-j) < sqrradius) and (k< sqrdim*lambda) then
           begin
             x[ii,jj]:= 1;
             k := k+1  ;
           end;
        end;
      end;
    else if showwarnings then ShowMessage('Randomize Patches nicht richtig aufgerufen');
  end;
end;


procedure outfile_initialize;
// initializes outfile
var parameters : string;
var time : string;
begin
  time := FormatDateTime('yy-mm-dd_hh-mm-ss_',Now);
  parameters := 'statistics\'+ time + statisticsfilename +'.ini';
  write_ini(parameters);
  AssignFile(outfile, 'statistics\'+ time + statisticsfilename + '.dat');
  rewrite(outfile);
  append(outfile);
  writeln(outfile, '# Remarks:  ' + remarks);
  writeln(outfile, '#');
  writeln(outfile, '#/ / MODEL PARAMETERS');
  writeln(outfile, '#/ / random seed =', RandSeed);
  writeln(outfile, '#/ / dim=',dim);
  writeln(outfile, '#/ / lambda=', lambda);
  writeln(outfile, '#/ / sigma=', sigma);
  writeln(outfile, '#/ / w=',w);
  writeln(outfile, '#/ / interactionlength=',interactionlength);
  writeln(outfile, '#/ / pricestability=',pricestability  );
  writeln(outfile, '#/ / desiredbenefit =', desiredbenefit);
  writeln(outfile, '#/ / Ecological Model:');
  writeln(outfile, '#/ / mortality = ' , mortality  );
  writeln(outfile, '#/ / rho = ' , rho );
  writeln(outfile, '#/ / alpha = ' , alpha  );
  writeln(outfile, '#/ / dispersalrate = ' , dispersalrate  );
  writeln(outfile, '#');
  // RUN
  writeln(outfile, '#/ / RUN PARAMETERS');
  writeln(outfile, '#/ / run settings');
  writeln(outfile, '#/ / trades = ', trades );
  writeln(outfile, '#/ / intensity = ', intensity );
  writeln(outfile, '#/ / maxt =' , maxt);
  writeln(outfile, '#/ / runs = ', runs);
  writeln(outfile, '#/ / serial = ', serial);
  writeln(outfile, '#/ / altruism = ', altruism);
  writeln(outfile, '#/ / altruismon = ', altruismon);
  writeln(outfile, '#/ / ecologicalmodelon = ', ecologicalmodelon);
  writeln(outfile, '#/ / economicmodelon = ', economicmodelon);
  writeln(outfile, '#');
  // VARIATION
  writeln(outfile, '#/ / VARIATION SETTINGS');
  writeln(outfile, '#/ / 1d');
  writeln(outfile, '#/ / varyaround =' , varyaround);
  writeln(outfile, '#/ / Steps =' , varysteps);
  writeln(outfile, '#/ / Range = ' , varyrange);
  writeln(outfile, '#/ / tries = ' , varytries);
  writeln(outfile, '#/ / 2d/ / ');
  writeln(outfile, '#/ / varyaround2d1 =' , varyaround2d1);
  writeln(outfile, '#/ / Steps2d1 =' , varysteps2d1);
  writeln(outfile, '#/ / Range2d1 = ' , varyrange2d1);
  writeln(outfile, '#/ / varyaround2d2 =' , varyaround2d2);
  writeln(outfile, '#/ / Steps2d2 =' , varysteps2d2);
  writeln(outfile, '#/ / Range2d2 = ' , varyrange2d2);
  writeln(outfile, '#/ / showvalues = ' , showvalues);
  writeln(outfile, '#/ / showwarnings = ' , showwarnings);
  writeln(outfile, '#/ / random seed = ' , randseed);
  writeln(outfile, '#');
  //writeln(outfile, 'loop #' , chr(9), 'xtot' ) ;
  //writeln(outfile);

  writeln(outfile, 'dimensions' , chr(9),
                   'lambda' , chr(9),
                   'w', chr(9),
                   'interactionlength', chr(9),
                   'neighbornumber', chr(9),
                   'sigma', chr(9),
                   'pricestability',chr(9),
                   'habitats', chr(9),'sdhabitat',chr(9),
                   'turnover', chr(9),'sdturnover',chr(9),
                   'connectivity', chr(9),'sdconnectivity',chr(9),
                   'totalcosts', chr(9),'sdtotalcosts',chr(9),
                   'avpayments',chr(9),'sdavpayments', chr(9),
                   'price',chr(9),'sdprice', chr(9),
                   'populationfrequency',chr(9),'sdpopulationfrequency', chr(9),
                   'survival after' + inttostr(maxt),chr(9),'empty'
                   );
  writeln(outfile, '#');
end;

procedure outfile_write ;
// writes outfile
begin
append(outfile);
writeln(outfile ,dim, chr(9), lambda, chr(9),w, chr(9), interactionlength , chr(9),neighbornumber , chr(9),sigma, chr(9),
pricestability,chr(9),
summary[4,1], chr(9),  summary[4,2] ,chr(9),
summary[1,1] , chr(9), summary[1,2] , chr(9),
summary[3,1],  chr(9) ,summary[3,2],  chr(9),
summary[2,1], chr(9),  summary[2,2], chr(9),
summary[5,1]*desiredbenefit ,chr(9), summary[5,2]*desiredbenefit  ,chr(9),
summary[5,1],chr(9),   summary[5,2],chr(9),
summary[6,1] , chr(9), summary[6,2] , chr(9),
summary[7,1] , chr(9), summary[7,2] );
end;


procedure costfield_write;
// writes the current cost distribution
var i,j : integer;
costfield : TextFile;
begin
  AssignFile(costfield, 'statistics\costfield.dat');
  rewrite(costfield);
  append(costfield);
  writeln(costfield, '# Remarks: Costfield' );
  for i:= 1 to dim do
  begin
        for j := 1 to dim do
        begin
          write(costfield, c[i,j]);
          if j < dim then write(costfield, chr(9));
        end;
  writeln(costfield, '');
  end;
  CloseFile(costfield);
end;

procedure costseries_write;
// writes the current cost distribution
var i,j : integer;
costfield : TextFile;
//costarray : array of double ;
begin
  ecologicalmodelon := false;
  AssignFile(costfield, 'statistics\costseries.dat');
  rewrite(costfield);
  append(costfield);
  writeln(costfield, '# Remarks: Costfield' );
  //setlength(costarray, maxt);
  for i:= 1 to runs do
  begin
        run(form1.RadioGroupCosts.Itemindex);
        writeln(costfield, c[15,15]);
  end;
  CloseFile(costfield);
end;





procedure logwindow_initialize;
// initializes the log window
begin
 if showvalues then
 begin
 form1.RunStatistics.TabWidth := 28;
     form1.RunStatistics.Items.insert(0,
     '#'+ chr(9) + 'time'+ #9 +
     '# of Runs= ' + #9 +
     'av.Run= '+ #9+
     'P= ' +    #9 +
     'av.P= ' +#9+
     'Prec.=' + #9 +
     'xtot=' + #9 +
     'mxtot=' + #9 +
     'Turnov=');
     form1.RunStatistics.Items.insert(1,' ') ;
 end;
end;

procedure logwindow_show (inp:integer);
// writes some standard outputs in the logwindow
begin
   if showvalues then
   begin
    form1.RunStatistics.Items.insert(0,
    inttostr(inp)+#9+FormatDateTime('hh:nn:ss:zzz', ( Now-storetime)) + 'ms'+ #9
    +inttostr(loops)+ #9 +
    floattostrF(loopsum/maxt/trades,ffFixed, 3,3)+ #9+
    floattostrF(price,ffFixed, 3,3)+ #9+
    floattostrF(mprice,ffFixed, 3,3)+#9+
    floattostrF(ztot/desiredbenefit-1,ffexponent,3,2)+ #9 +
    inttostr(xtot)+ #9 +
    floattostrF(mxtot, ffFixed,4,2)+ #9 +
    floattostrF(turnoverrate,fffixed,3,2)
    );
   end;
end;

procedure picture_draw(inp,temp:integer);
// shows the graphics
var i,j,l,k: integer;
begin
  if showvalues then
  begin
    bild:=TBitmap.Create;
    try
    bild.width:=zoom*dim;
    bild.height:=zoom*dim;

    case temp of
    0:
    begin
    for i:=1 to dim do for j:=1 to dim do
      // habitats
      if (deathmatrix[i,j]= 1) and form1.checkboxshowdeaths.checked then for k:=1 to zoom do for l:=1 to zoom do bild.canvas.pixels[zoom*(i-1)+k,zoom*(j-1)+l]:=clred else
      if x[i,j]=1 then
      begin
        if (p[i,j]=1) and ecologicalmodelon then
          for k:=1 to zoom do for l:=1 to zoom do bild.canvas.pixels[zoom*(i-1)+k,zoom*(j-1)+l]:=clgreen
        else if (p[i,j]=2) and ecologicalmodelon then
          for k:=1 to zoom do for l:=1 to zoom do bild.canvas.pixels[zoom*(i-1)+k,zoom*(j-1)+l]:=clred
        else if form1.checkboxshownewlycreatedhabitats.checked and (xold[i,j]=0) then
          for k:=1 to zoom do for l:=1 to zoom do bild.canvas.pixels[zoom*(i-1)+k,zoom*(j-1)+l]:=cldkgray
        else for k:=1 to zoom do for l:=1 to zoom do bild.canvas.pixels[zoom*(i-1)+k,zoom*(j-1)+l]:=clblack;
      end
      else if x[i,j]=0 then
      begin
        if (xold[i,j]=1) and form1.checkboxshowdestroyedhabitats.checked then
          for k:=1 to zoom do for l:=1 to zoom do bild.canvas.pixels[zoom*(i-1)+k,zoom*(j-1)+l]:=clred
        else for k:=1 to zoom do for l:=1 to zoom do bild.canvas.pixels[zoom*(i-1)+k,zoom*(j-1)+l]:=clwhite;
      end
      else if showwarnings then showmessage('Fehler in picture_draw');
    end;
    1:
    begin
      for i:=1 to dim do for j:=1 to dim do
          if xtemp[i,j]=1 then for k:=1 to zoom do for l:=1 to zoom do bild.canvas.pixels[zoom*(i-1)+k,zoom*(j-1)+l]:=clblack
          else for k:=1 to zoom do for l:=1 to zoom do bild.canvas.pixels[zoom*(i-1)+k,zoom*(j-1)+l]:=clwhite;
    end;
    end;
    form1.image1.picture.assign(bild);
    if (inp > 0) and form1.checkboxsavevideo.checked then
    begin
      forcedirectories(ExtractFilePath(ParamStr(0))+'\bilder');
      bild.SaveToFile(ExtractFilePath(ParamStr(0))+'\bilder\'+inttostr(inp)+'.bmp');
    end;
    form1.refresh;
    finally
    bild.free;
    end;
  end;  
end;

procedure parameters_show;
// shows the present values of the parameters in the forms
begin
  if showvalues then
  begin
    enableforms;
    form1.edit1.text:=floattostr(w);
    form1.EditLength.text:= floattostr(interactionlength);
    form1.LabelNeighbors.caption := inttostr(neighbornumber) + ' neighbors';
    form1.edit2.text:=floattostr(sigma);
    form1.edit3.text:=inttostr(maxt);
    form1.edit7.text:=inttostr(dim);
    form1.edit9.text:=floattostr(lambda);
    form1.edit10.text:=inttostr(runs);
    form1.edit11.text:=floattostr(pricestability);
    form1.varound.text:=floattostr(varyaround);
    form1.vrange.text:=floattostr(varyrange);
    form1.vsteps.text:=floattostr(varysteps);
    form1.traderuns.text:=inttostr(trades);
    form1.variationtries.text:=inttostr(varytries);
    form1.varound2d1.text:=floattostr(varyaround2d1);
    form1.vrange2d1.text:=floattostr(varyrange2d1);
    form1.vsteps2d1.text:=floattostr(varysteps2d1);
    form1.varound2d2.text:=floattostr(varyaround2d2);
    form1.vrange2d2.text:=floattostr(varyrange2d2);
    form1.vsteps2d2.text:=floattostr(varysteps2d2);
    form1.editintensity.text := floattostr(intensity);
    form1.editprecision.text := floattostr(precision);
    form1.Editmortality.text := floattostr(mortality) ;
    form1.Editmortalitycorrelation.text := floattostr(rho) ;
    form1.Editdispersalrate.text := floattostr(dispersalrate) ;
    form1.Editdispersaldistance.text := floattostr(alpha) ;
    form1.editaltruism.Text := floattostr(altruism);
    form1.CheckBoxAltruismon.Checked := altruismon;
    form1.filename.text:= statisticsfilename;
    form1.EditInifilename.text := extractfilename(inifilename) ;
    Form1.EditOptimizationPercentage.text := floattostr(OptimizationPercentage);
    Form1.CheckBoxoptimizationon.checked := optimizationon;
    form1.RadioGroupVariationChoice.itemindex := variationchoice;
    form1.EditReduce.text:= inttostr(reduce);

    form1.CheckBoxEcologicalmodelon.Checked := ecologicalmodelon;
    form1.CheckBoxEconomicmodelon.Checked := economicmodelon;
    form1.CheckBoxShowWarnings.Checked := showwarnings;
    form1.editzoom.text:=inttostr(zoom);
    form1.radiogroupserial.itemindex := serial;
    form1.RadioGroupCosts.itemindex:= costvariation;
    form1.RadioGroupStartDistribution.itemindex:=startdistribution;
    form1.RadiogroupTradeConstraints.ItemIndex := tradeconstraints;
    form1.EditCostSize.text := inttostr(costsize);
    form1.CheckBoxIndependentrunon.checked := independentrunon;

    form1.EditGTcc.text := floattostr(cc);
    form1.EditGTdc.text := floattostr(dc);
    form1.EditGTcd.text := floattostr(cd);
    form1.EditGTdd.text := floattostr(dd);

    form1.Editvradius.text := floattostr(vradius);
    form1.Editmradius.text := floattostr(mradius);
    form1.Editdispersalpreference.text := floattostr(dispersalpreference);

    //optimization
    form1.editTemperatureDecay.text := floattostr(temperaturedecay);
    form1.editstarttemperature.text := floattostr(starttemperature);
  end;
end;

procedure values_show;
// shows the present calculated statistics
var i : integer;
begin
  if showvalues then
  begin
    case tradeconstraints of
    0:
      begin
        form1.edit4.text:=floattostrF(desiredBenefit, fffixed,4,3);
        form1.goaldiff.text:=floattostrF((ztot/desiredbenefit-1)*100, fffixed,4,3)+'%';
      end;
    1:
      begin
        form1.goaldiff.text:=floattostrF( (totalcosts/lambda-1)*100, fffixed,4,3)+'%';
        form1.edit4.text:=floattostrF(lambda, fffixed,4,3);
      end;
    end;
    form1.edit5.text:=floattostrF(xtot/sqrdim*100, fffixed,4,3)+'%';
    form1.edit6.text:=floattostrF(connectivity*100, fffixed,4,3) +'%';
    form1.edit8.text:=floattostrF(price, fffixed,4,3);
    form1.turnover.text:=floattostrF(turnoverrate, fffixed,4,3);
    form1.costs.text:=floattostrF(totalcosts, fffixed,4,3);

    form1.Neighboursshow.Items.Clear;
    form1.Neighboursshow.Items.Add('Nachbarn Anteil');
    form1.Neighboursshow.Items.Add('------------------');
    form1.Series5.Clear;
    for  i:=0 to neighbornumber do
    begin
      form1.series5.AddBar(neighbourstatistics[i],'', clBlue );
      form1.Neighboursshow.Items.Add(inttostr(i) + ' -- ' + floattostrF(neighbourstatistics [i],fffixed,4,3)+'%' );
    end;
    form1.Neighboursshow.Items.Add('av Neighbours');
    form1.Neighboursshow.Items.Add('------------------');
    form1.Neighboursshow.Items.Add( floattostrF(averageneigbours, fffixed,4,3) );
    form1.costlevel.text:=(floattostrF(averagecostlevel, fffixed,4,3));
    form1.editaveragepayments.text := (floattostrF(price*desiredbenefit, fffixed,4,3));
    form1.ShowPopulatedHabitats.text := floattostrF(numberofpopulatedhabitats/xtot*100, fffixed,4,3)+'%';
    form1.ShowColonizationProbability.text := floattostr(mcolonizationsum);
    form1.showsdcosts.text := (floattostrF(sdcosts, fffixed,4,3));
    form1.showmcostchange.text := (floattostrF(mcostchange, fffixed,4,3));

    Form1.Memo1.Clear;
    Form1.Memo1.text := remarks;
  end;  
end;

procedure read_ini;
// reads the ini file
var ini: TIniFile;
begin
  ini:=TIniFile.Create(inifilename);
  try
        // parameters
  
        lambda:= ini.ReadFloat('Parameter','lambda',0.15);
        dim:= ini.ReadInteger('Parameter','dim',50);
        w:=ini.ReadFloat('Parameter','w',0.32);
        interactionlength:= ini.ReadFloat('Parameter','interactionlength',1);
        sigma:=ini.ReadFloat('Parameter','sigma',0.6);
        pricestability:= ini.ReadFloat('Parameter','pricestability',0.01);
        altruism:= ini.ReadFloat('Parameter','altruism',0.5);


        // Ecological Model

        rho:=ini.ReadFloat('Ecological Model','rho',2);
        dispersalrate:=ini.ReadFloat('Ecological Model','dispersalrate',0.2);
        alpha:=ini.ReadFloat('Ecological Model','alpha',3);
        mortality:=ini.ReadFloat('Ecological Model','mortaltiy',0.2);
        vradius:=ini.ReadFloat('Ecological Model','vradius',0);
        mradius:=ini.ReadFloat('Ecological Model','mradius',2);

        //run

        maxt:=ini.ReadInteger('Run','maxt',200);
        runs:= ini.ReadInteger('Run','runs',5);
        trades := ini.ReadInteger('Run','trades',1);
        price :=ini.ReadFloat('Run','price',0.5);
        intensity := ini.ReadInteger('Run','intensity',1);
        precision :=ini.ReadFloat('Run','precision',0.01);
        serial :=ini.Readinteger('Run','serial',0);
        costvariation := ini.Readinteger('Run','costvariation',0);
        startdistribution :=ini.Readinteger('Run','startdistribution',0);

        // variation

        varyRange :=ini.ReadFloat('Variation','varyrange',0.15);
        varysteps := ini.ReadInteger('Variation','varysteps',5);
        varytries := ini.ReadInteger('Variation','varytries',5);
        varyaround :=ini.ReadFloat('Variation','varyaround',0.5);

        varyaround2d1 :=ini.ReadFloat('Variation','varyaround2d1',0.5);
        varyRange2d1 :=ini.ReadFloat('Variation','varyRange2d1',0.15);
        varysteps2d1 := ini.ReadInteger('Variation','varysteps2d1',5);

        varyaround2d2 :=ini.ReadFloat('Variation','varyaround2d2',0.5);
        varyRange2d2 :=ini.ReadFloat('Variation','varyRange2d2',0.15);
        varysteps2d2 := ini.ReadInteger('Variation','varysteps2d2',5);

        variationparameterchoice := ini.ReadInteger('Variation','variationparameterchoice',0);
        variationchoice := ini.ReadInteger('Variation','variationchoice',0);
        statisticsfilename := ini.ReadString('Variation','statisticsfilename','statistics.dat');
        remarks := ini.ReadString('Variation','remarks','no remarks');
        reduce := ini.ReadInteger('Variation','reduce',1);
        // options

        showvalues := ini.ReadBool('Options','showvalues',false);
        ecologicalmodelon := ini.ReadBool('Options','ecologicalmodelon',false);
        economicmodelon := ini.ReadBool('Options','economicmodelon',true);
        showwarnings := ini.ReadBool('Options','showwarnings',true);
        zoom := ini.ReadInteger('Options','zoom',3);
        altruismon := ini.ReadBool('Options','altruismon',false);
        tradeconstraints := ini.readinteger('Options','tradeconstraints',0);
        independentrunon := ini.ReadBool('Options','independentrunon',false);
        costsize := ini.readinteger('Options','costsize',1);

        // optimization

        optimizationon := ini.ReadBool('Optimization','OptimizationOn',false);
        OptimizationPercentage :=ini.ReadFloat('Optimization','OptimizationPercentage',0.1);
        StartTemperature :=ini.ReadFloat('Optimization','StartTemperature',30);
        TemperatureDecay :=ini.ReadFloat('Optimization','TemperatureDecay',0.01);

        //game theory

        cc :=ini.ReadFloat('GameTheory','cc',0.8);
        cd :=ini.ReadFloat('GameTheory','cd',0.1);
        dc :=ini.ReadFloat('GameTheory','dc',1);
        dd :=ini.ReadFloat('GameTheory','dd',0.2);
  finally
    ini.free;
  end;
  constants_calculate;
end;

procedure write_ini(filename:string);
// writes the ini file
var ini: TIniFile;
begin
  constants_calculate;
  ini:=TIniFile.Create(filename);
  try
        // economic parameters

        ini.writeFloat('Parameter','lambda',lambda);
        ini.writeInteger('Parameter','dim',dim);
        ini.writeFloat('Parameter','w',w);
        ini.writeFloat('Parameter','interactionlength',interactionlength);
        ini.writeFloat('Parameter','sigma',sigma);
        ini.writeFloat('Parameter','pricestability',pricestability);
        ini.writeFloat('Parameter','altruism',altruism);

        // Ecological Model

        ini.writeFloat('Ecological Model','mortaltiy',mortality);
        ini.writeFloat('Ecological Model','rho',rho);
        ini.writeFloat('Ecological Model','dispersalrate',dispersalrate);
        ini.writeFloat('Ecological Model','alpha',alpha);
        ini.writeFloat('Ecological Model','vradius',vradius);
        ini.writeFloat('Ecological Model','mradius',mradius);

        //run

        ini.writeInteger('Run','maxt',maxt);
        ini.writeInteger('Run','runs',runs);
        ini.writeInteger('Run','trades',trades);
        ini.writeFloat('Run','price',price);
        ini.writeInteger('Run','intensity',intensity );
        ini.ReadBool('Run','showwarnings',showwarnings);
        ini.writeFloat('Run','precision',precision);
        ini.writeInteger('Run','serial',serial);
        ini.writeinteger('Run','costvariation',costvariation);
        ini.writeinteger('Run','startdistribution',startdistribution);

        // variation

        ini.writeFloat('Variation','varyrange',varyRange);
        ini.writeInteger('Variation','varysteps', varysteps);
        ini.writeFloat('Variation','varyaround',varyaround);

        ini.writeFloat('Variation','varyaround2d1',varyaround2d1);
        ini.writeFloat('Variation','varyRange2d1',varyRange2d1);
        ini.writeInteger('Variation','varysteps2d1',varysteps2d1 );

        ini.writeFloat('Variation','varyaround2d2',varyaround2d2);
        ini.writeFloat('Variation','varyRange2d2',varyRange2d2);
        ini.writeInteger('Variation','varysteps2d2', varysteps2d2);

        ini.writeInteger('Variation','variationparameterchoice',variationparameterchoice);
        ini.writeInteger('Variation','variationchoice',variationchoice);
        ini.writeInteger('Variation','varytries',varytries);
        ini.writeString('Variation','statisticsfilename',statisticsfilename);
        ini.writeString('Variation','remarks',remarks);
        ini.writeInteger('Variation','reduce',reduce);
        //options

        ini.writeBool('Options','showvalues',showvalues);
        ini.writeBool('Options','ecologicalModelOn',ecologicalmodelon);
        ini.writeBool('Options','economicmodelon',economicmodelon);
        ini.writeBool('Options','showwarnings',showwarnings);
        ini.writeInteger('Options','zoom',zoom);
        ini.writebool('Options','altruismon',altruismon) ;
        ini.writeInteger('Options','tradeconstraints',tradeconstraints);
        ini.writebool('Options','independentrunon',independentrunon) ;
        ini.writeInteger('Options','costsize',costsize);

        // optimization

        ini.writebool('Optimization','OptimizationOn',optimizationon) ;
        ini.writeFloat('Optimization','OptimizationPercentage',OptimizationPercentage);
        ini.writeFloat('Optimization','StartTemperature',StartTemperature);
        ini.writeFloat('Optimization','TemperatureDecay',TemperatureDecay);

        // Game Theory

        ini.WriteFloat('GameTheory','cc',cc);
        ini.WriteFloat('GameTheory','cd',cd);
        ini.WriteFloat('GameTheory','dc',dc);
        ini.WriteFloat('GameTheory','dd',dd);
  finally
    ini.free;
  end;
end;

function compareProfit(item1 : pointer; item2 : pointer) : Integer;
// used for the optimization
// implements a sorting function for list
// item gets a smaller value if ecological benefit / costs is smaller
// = worst items are first
var
  x1,x2,y1,y2:integer;
  coordinate1,coordinate2 : coordinatepointer ;
begin
  coordinate1 := item1;
  coordinate2 := item2;

  x1 := coordinate1^[0];
  x2 := coordinate2^[0];
  y1 := coordinate1^[1];
  y2 := coordinate2^[1];

  if   ( ecological_benefit_alt(mtemp[x1,y1])/c[x1,y1])>
  ( ecological_benefit_alt(mtemp[x2,y2])/c[x2,y2])
  then Result := -1
  else if ( (ecological_benefit_alt(mtemp[x1,y1])/c[x1,y1])<
  ( ecological_benefit_alt(mtemp[x2,y2])/c[x2,y2] ) )
  then Result := 1
  else Result := 0;
end;

function transitionProbability(current, flip, temperature: double): double ;
// optimization function
// calculates the transition probability to flip at a given temperature and
// the two energies of the present and the flipped state;
begin
   case optimizationAlgorithm of
     0:
     begin
       if (flip - current) > 0 then result := exp(- (flip-current)/temperature) else result := 1;
       //if (flip - current) > 0 then result:= 0 else result := 1;
     end;
   end;
end;


function flip(inp : integer) : integer;
// flips 1 to 0 and 0 to 1
begin
if inp = 1 then result := 0 else
if inp = 0 then result := 1 ;
end;

procedure optimalTrade;
// optimization main function
var i,j,ii,jj,k,l,kk,ll: integer;

    uplist, downlist:TList;
    Koordinate:^coordinate;

    // SA variables
    benefitcurrent,benefitflip, benefitup, benefitdown: double;
    temperature : double;  // current temperature

    // to store the best solution found
    rand,tran : double;
    ctot: double ;
    xOpt : intmatrix;
    cbopt : double;
    
begin

  // initializations

  temperature := startTemperature;
  cbopt := 1000000;  // set arbitrarily high so that the first solution is accepted
  for i:=1 to dim do  for j:=1 to dim do xtemp[i,j]:=x[i,j]; // copy mtemp
  calculateNumberOfNeighbors_xtemp_mtemp;

  if optimizationshow then
  begin
   form2.series1.clear;
   form2.series2.clear;
   form2.series3.clear;
  end;

  // SA Procedure

  for k:=1 to trades do
  begin

       // decision process  ---------------------


        xtot:=0;
        ztot:=0;
        ctot  :=0;
        for i:=1 to dim do  for j:=1 to dim do
        begin
            xtot :=xtot+xtemp[i,j];
            ztot:=ztot+xtemp[i,j]*(ecological_benefit(mtemp[i,j]));
            ctot :=ctot+xtemp[i,j]*c[i,j];
        end;
     {
        for l:= 0 to ceil((sqrdim)*optimizationpercentage) do // takes a percentage of random grid cells without repetition
        begin
          i := random(dim)+1; j := random(dim)+1;
          benefitcurrent := ctot/ztot;
          benefitflip := c[i,j]/(marginal_ecological_benefit(mtemp[i,j]));
          rand := random;
          tran := transitionProbability(benefitcurrent, benefitflip, temperature);
          if rand < tran then xtemp[i,j]:=1 else xtemp[i,j]:= 0;
           calculatesinglemtemp(i,j);
       end;    }

       for l:= 0 to ceil((sqrdim)*optimizationpercentage) do // takes a percentage of random grid cells without repetition
        begin
          i := random(dim)+1; j := random(dim)+1;
          benefitcurrent := ctot/ztot;
          benefitflip := (flip(x[i,j])*c[i,j]/marginal_ecological_benefit(mtemp[i,j])-x[i,j]*c[i,j]/marginal_ecological_benefit(mtemp[i,j])+sqrdim*benefitcurrent)/sqrdim;
          rand := random;
          tran := transitionProbability(benefitcurrent, benefitflip, temperature);
          if rand < tran then xtemp[i,j]:=1 else xtemp[i,j]:= 0;
           calculatesinglemtemp(i,j);
       end;


        calculateNumberOfNeighbors_xtemp_mtemp;

        // check that constraints are met
        xtot:=0;
        ztot:=0;
        ctot  :=0;
        for i:=1 to dim do  for j:=1 to dim do
        begin
            xtot :=xtot+xtemp[i,j];
        end;


        uplist:=Tlist.Create;
        downlist:=Tlist.Create;
        uplist.capacity:=xtot+1;
        downlist.capacity:=sqrdim-xtot+1;
        try
          for i:=1 to dim do  for j:=1 to dim do
          begin
           New(Koordinate);
           Koordinate^[0]:= i; Koordinate^[1]:= j ;
           if xtemp[i,j]= 1 then
           uplist.add(Koordinate) else downlist.add(Koordinate);
          end;
          uplist.sort(compareProfit);
          downlist.sort(compareProfit);

          ztot:=0;
          calculateNumberOfNeighbors_xtemp_mtemp;
          for i:=1 to dim do  for j:=1 to dim do
              begin
                ztot:=ztot+xtemp[i,j]*(ecological_benefit(mtemp[i,j]));
              end;

              if  (ztot-desiredbenefit*sqrdim) > 0 then
              begin
                l:=(uplist.count-1);
                repeat
                  ztot:=0;
                  for i:=1 to dim do  for j:=1 to dim do
                  begin
                    ztot:=ztot+xtemp[i,j]*(ecological_benefit(mtemp[i,j]));
                  end;
                  if ztot-desiredbenefit*sqrdim < 0 then break else
                  Koordinate:= uplist.Items[l] ;
                  i:= Koordinate^[0]; j:= Koordinate^[1];
                  xtemp[i,j]:= 0;
                  calculatesinglemtemp(i,j);
                  l:= l -1 ;
                until false ;
              end else
              if ( ztot-desiredbenefit*sqrdim) <0 then
              begin
                l:=0;
                repeat
                  ztot:=0;
                  for i:=1 to dim do  for j:=1 to dim do
                  begin
                    ztot:=ztot+xtemp[i,j]*(ecological_benefit(mtemp[i,j]));
                  end;
                  if ztot-desiredbenefit*sqrdim > 0 then break else
                  Koordinate:= downlist.Items[l] ;
                  i:= Koordinate^[0]; j:= Koordinate^[1];
                  xtemp[i,j]:= 1;
                  calculatesinglemtemp(i,j);
                  l:= l +1 ;
                until false;
              end;
        finally
          for l := (uplist.count-1) downto 0 do
          begin
           Koordinate:= uplist.Items[l];
           Dispose(Koordinate);

          end;
          for l := (downlist.count-1) downto 0 do
          begin
           Koordinate:= downlist.Items[l];
           Dispose(Koordinate);
          end;
          while uplist.Count>0 do Uplist.Delete(0);
          while downlist.Count>0 do downlist.Delete(0);
          uplist.Free;
          downlist.Free;
        end;

        //store the best solution
        xtot:=0;
        ztot:=0;
        ctot  :=0;
        for i:=1 to dim do  for j:=1 to dim do
          begin
            xtot :=xtot+xtemp[i,j];
            ztot:=ztot+xtemp[i,j]*(ecological_benefit(mtemp[i,j]));
            ctot :=ctot+xtemp[i,j]*c[i,j];
          end;
        if ctot/ztot < cbopt then
        begin
           cbopt := ctot/ztot;
           for i:=1 to dim do  for j:=1 to dim do
           begin
            xOpt[i,j]:= xtemp[i,j];
            x[i,j]:= xtemp[i,j];
           end;
        end;

       //update temperature

       temperature := temperature*(1-temperaturedecay);
       if optimizationshow and (k mod 5 = 0) and (k>(trades/8)) then
       begin
       picture_draw(0,1);
       calculations_do;
       values_show;
       Form2.Show;
       with form2.Series2 do AddXY(k,temperature);
       with form2.Series1 do AddXY(k,totalcosts);
       with form2.Series3 do AddXY(k,cbopt);
       form2.refresh;
       end;
  end; // end trades loop

  // write the best solution to x
  for i:=1 to dim do  for j:=1 to dim do x[i,j]:=xopt[i,j];

  if optimizationshow then
  begin
  //picture_draw(0,0);
  end;
end;


procedure trade;
// implements the trading of the market
// serial trade (random people trade, after a while the knowledge is updated)

var i,j,ii,jj,k,l,n,ran,loop: integer;
    price_factor : double ;
    ztot_old : double;
    selection: array of integer;
    partition: array of integer;
    var Liste:TList;
    Zahl:^Integer;
    pre : double;
    neighborsij : neighborcoordinates ;
    benefitmit,benefitohne : double;

begin
   // initializations --------------------------------------
   // random without repetition init
   if (serial = 1) or (serial = 2) then
   begin
     SetLength(selection, sqrdim);
     SetLength(partition, intensity);
     Liste:=TList.Create;
     Liste.capacity:=sqrdim;
     for i := 0 to sqrdim-1 do
     begin
       New(Zahl);
       Zahl^:=i;
       Liste.Add(Zahl);
     end;
     for i := 0 to sqrdim-1 do
     begin
       ran := random(sqrdim-i);
       Zahl:=Liste.Items[ran];
       selection[i] := Zahl^  ;
       Dispose(Zahl);
       Liste.Delete(ran);
     end;
     if Liste.Count <> 0 then ShowMessage('Problem in tradeserial');
     Liste.Free;
   end;

   // partition init
   if serial = 1 then
   begin
     i:= intensity -1 ;
     j:= round (sqrdim/intensity);
     while i > 0 do
     begin
     partition[i] := j  ;
     i := i-1;
     end;
     partition[0] := sqrdim - j * (intensity -1);
   end;

   // rest init
   calculateNumberOfNeighbors_x_m;
   price_factor := price ;
   loops :=0 ; ztot:=0;
   case tradeconstraints of
   0:
   begin
     for i:=1 to dim do  for j:=1 to dim do  // initialize ztot and help matrices
     begin
     ztot:=ztot+x[i,j]*(ecological_benefit(m[i,j]));
     end;
     ztot := ztot/sqrdim;
   end;
   1:
   begin
     ztot := totalcosts;
   end;
   else if showwarnings then showmessage('Fehler in Run Tradeconstraints')
   end;

   // decision process  ---------------------

   pre := precision;
   repeat
      case serial of
      0:         // Paralled Decisions ... CAREFULL: XTEMP IS NOT INITIALIZED!!!
        begin
          for i:=1 to dim do  for j:=1 to dim do xtemp[i,j]:=x[i,j];
          calculateNumberOfNeighbors_xtemp_mtemp;
          for i:=1 to dim do  for j:=1 to dim do
          begin
            if altruismon then
            begin
              if -c[i,j]+price_factor*(ecological_benefit_alt(m[i,j]))>0 then
              xtemp[i,j]:=1
              else xtemp[i,j]:=0;
            end
            else
            begin
              if -c[i,j]+price_factor*(ecological_benefit(m[i,j]))>0 then
              xtemp[i,j]:=1
              else xtemp[i,j]:=0;
            end;
          end;
        end;
      1:        // partly parallel, partly sequential decisions
        begin
          loop:=0 ;   // sets loop counter to 0 ; Loop counts the runs inbetween 2 loops
          for i:=1 to dim do  for j:=1 to dim do xtemp[i,j]:=x[i,j]; // so that the right mtemp is initialized
          for l:= 0 to (intensity -1) do
          begin
            calculateNumberOfNeighbors_xtemp_mtemp;
            for k:=1 to partition[l] do
            begin
              i := (selection[loop] div dim ) +1 ; j:= (selection[loop] mod dim )+ 1 ;
              if altruismon then
              begin
                if -c[i,j]+price_factor*(ecological_benefit_alt(m[i,j]))>0 then
                xtemp[i,j]:=1
                else xtemp[i,j]:=0;
              end
              else
              begin
                if -c[i,j]+price_factor*(ecological_benefit(mtemp[i,j]))>0 then
                xtemp[i,j]:=1
                else xtemp[i,j]:=0;
              end;
              loop := loop +1 ;
            end;
          end;
        end;
      2:              // sequential
        begin
          for i:=1 to dim do  for j:=1 to dim do xtemp[i,j]:=x[i,j]; // copy xtemp
          calculateNumberOfNeighbors_xtemp_mtemp;
          for l:= 0 to sqrdim-1 do
          begin
            i := (selection[l] div dim ) +1 ; j:= (selection[l] mod dim )+ 1 ;
            if altruismon then
            begin
              if -c[i,j]+price_factor*(ecological_benefit_alt(m[i,j]))>0 then
              xtemp[i,j]:=1
              else xtemp[i,j]:=0;
            end
            else
            begin
              if -c[i,j]+price_factor*(ecological_benefit(mtemp[i,j]))>0 then
              xtemp[i,j]:=1
              else xtemp[i,j]:=0;
            end;
            //calculateNumberOfNeighbors_xtemp_mtemp;    // as an alternative for the singlemtemp function
            if xtemp[i,j] <> x[i,j] then calculatesinglemtemp(i,j);
          end;
        end;
        3:         // follow costs only 
        begin
          for i:=1 to dim do  for j:=1 to dim do
          begin
              if -c[i,j]+price_factor > 0 then
              xtemp[i,j]:=1
              else xtemp[i,j]:=0;

          end;
        end;
      else if showwarnings then showmessage('Wrong serial in Tradeserial');
      end;

      // Check if goal is hit --------------------------

      case tradeconstraints of
      1:
        begin
          ztot_old := ztot; ztot:=0;
          for i:=1 to dim do  for j:=1 to dim do if xtemp[i,j]= 1 then ztot:=ztot+c[i,j];
          ztot := ztot/sqrdim;
          price := price_factor;
          if abs(ztot/lambda -1) < pre then break else
          price_factor:=price_factor - 0.03*sign(ztot/lambda -1)*random ;//+ 0.1*round(ztot/sqrdim*10) ;
          if abs(ztot/lambda-1) > abs(ztot_old/lambda-1) then price_factor:=price;
          loops:= loops +1 ;
          if loops > 500 then
          begin
          pre := pre * 2 ;
          if showwarnings then form1.RunStatistics.Items.insert(0, 'Precision wurde auf ' + floattostr(pre) + 'angepasst');
          loops := 200;
          end;
        end;
      0:
        begin
          ztot_old := ztot; ztot:=0;
          calculateNumberOfNeighbors_xtemp_mtemp;
          for i:=1 to dim do  for j:=1 to dim do
          begin
            ztot:=ztot+xtemp[i,j]*(ecological_benefit(mtemp[i,j]));
          end;
          ztot := ztot/sqrdim;
          price := price_factor;
          if abs(ztot/desiredbenefit-1) < pre then break else
          price_factor:=price_factor - 0.03*sign(ztot/desiredbenefit-1)*random ;//+ 0.1*round(ztot/sqrdim*10) ;
          if abs(ztot) > abs(ztot_old) then price_factor:=price;
          loops:= loops +1 ;
          if loops > 250 then
          begin
          pre := pre * 2 ;
          if showwarnings then form1.RunStatistics.Items.insert(0, 'Precision wurde auf ' + floattostr(pre) + 'angepasst');
          loops := 200;
          end ;
        end;
      else if showwarnings then showmessage('Fehler in Run Tradeconstraints')
      end;
   until False;
   for i:=1 to dim do  for j:=1 to dim do x[i,j]:=xtemp[i,j];
end;

// -----------------------------------------------------------------------------------
// Ecological Model

procedure initializePopulation;
// initializes the population of the ecological model
var i,j : integer;
    fractionofpopulatedhabitats: real;
    //minimumneighbours  : integer;
begin
  // variables
  fractionofpopulatedhabitats:=0.6;
  //minimumneighbours:= 3 ;
  for i:=1 to dim do  for j:=1 to dim do
  begin
  if (x[i,j]= 1) and ( random < fractionofpopulatedhabitats) {and (m[i,j]>= minimumneighbours)}
    then p[i,j]:= 1
    else p[i,j]:= 0
  end;
  ecocounter := 0;
end;


procedure initializeEcoStatistics(length : integer);
// sets lengts for statistics of the ecological model
begin
setlength(mortalitystatistics, length +1) ;
setlength(populationstatistics, length +1) ;
end;




procedure createPopulationList;
// lists of populated and unpopulated pathces
var i,j,l,n : integer;
begin
  l:=0; n:=0;
  for i:=1 to dim do  for j:=1 to dim do
  begin
    if p[i,j]=1 then
    begin
      populatedhabitats[l,0]:=i;
      populatedhabitats[l,1]:=j;
      l:=l+1;
    end
    else if (x[i,j]=1) and (p[i,j]= 0) then
    begin
      unpopulatedhabitats[n,0]:=i;
      unpopulatedhabitats[n,1]:=j;
      n:=n+1;
    end;
  end;
  numberofpopulatedhabitats:= l;
  numberofunpopulatedhabitats := n;
end;

procedure mortalitycircles(inp:integer);
// implements correlated environmental stochasticity
var death,pop: integer ;
    ii,jj,i,j,iii,jjj : integer;
    radius : double;
    sqrradius, sqrmradius : double;
    presentmortality:double;
    lognormalmu, lognormalsigma: double; // input for the gaussian function to produce lognormal distr.
begin
  death := 0 ;
  sqrmradius := sqr(mradius);
  //lognormalmu:= ln(mortality)- 0.5* ln(1+rho/sqr(mortality));
  //lognormalsigma := sqrt(ln(1+rho/sqr(mortality)));
  // presentmortality := exp(randg(lognormalmu, lognormalsigma)) ;

  presentmortality := mortality;

  //mortalitystatistics[ecocounter]:= presentmortality ;
  case inp of
    1:
    begin
      while death < (numberofpopulatedhabitats*presentmortality - (numberofpopulatedhabitats/xtot*sqrmradius*pi/2))  do
      begin
        lognormalmu:= ln(mradius)- 0.5* ln(1+vradius/sqrmradius);
        lognormalsigma := sqrt(ln(1+vradius/sqrmradius));
        radius := exp(randg(lognormalmu, lognormalsigma));
        sqrradius := sqr(radius);
        i := random(dim) +dim + 1 ; j := random(dim) +dim + 1;
        for ii:=i-round(radius) to i + round(radius) do for jj:=j-round(radius) to j + round(radius) do
        begin
           if sqr(ii-i) + sqr(jj-j) < sqrradius then
           begin
             iii := ii mod dim ; if iii = 0 then iii := dim;
             jjj := jj mod dim ; if jjj = 0 then jjj := dim;
             if p[iii, jjj] = 1 then
             begin
               p[iii, jjj]:= 0;
               death := death +1;
             end;
           end;
        end;
      end;
    end;
    2:
    begin
      while death < (sqrdim*presentmortality - (sqrmradius*pi/2)) do
      begin
        lognormalmu:= ln(mradius)- 0.5* ln(1+vradius/sqrmradius);
        lognormalsigma := sqrt(ln(1+vradius/sqrmradius));
        radius := exp(randg(lognormalmu, lognormalsigma));
        sqrradius := sqr(radius);
        i := random(dim) +dim + 1 ; j := random(dim) +dim + 1;
        for ii:=i-round(radius) to i + round(radius) do for jj:=j-round(radius) to j + round(radius) do
        begin
           if sqr(ii-i) + sqr(jj-j) < sqrradius then
           begin
             iii := ii mod dim ; if iii = 0 then iii := dim;
             jjj := jj mod dim ; if jjj = 0 then jjj := dim;
             if p[iii, jjj] <> 2 then
             begin
               p[iii, jjj]:= 2;
               death := death +1;
             end;
           end;
        end;
      end;


      pop := 0 ;
      for i:=1 to dim do  for j:=1 to dim do
      begin
        deathmatrix[i,j]:= 0;
        if p[i,j]= 1 then pop := pop +1;
        if p[i,j]= 2 then
        begin
          p[i,j]:= 0;
          deathmatrix[i,j]:= 1;
        end;
      end;
    death := numberofpopulatedhabitats-pop;
    end;
    else if showwarnings then showmessage('Fehler in mortalitycircles');
  end;
  if showvalues then
  begin
        if (numberofpopulatedhabitats > 0) then  mortalitystatistics[ecocounter]:= death/numberofpopulatedhabitats else mortalitystatistics[ecocounter]:= 0 ;
  end;
end;


function dispersalProbability (inp:double):double;
// imput is the euklidian distace between two patches
// dispersalrate controlls the amount of dispesal
// dispersal assumes a preference of 1 over sqr(distance) for close habitats
// dispersal is normalized with 1/alpha ... maybe this should be modified to a fraction of that to account for deaths during dispersal
// !!! The function needs additional normalization with sum 1 over d^2 to normalize the dispersal !!!
begin
  //if inp > 0 then result := (1/sqr(inp)) * (exp((-1)*(inp-1)/alpha))
  //else result := (exp((-1)*(inp-1)/alpha));
  result := exp(- (inp) / alpha)
end;

procedure disperse(i,j : integer; disprate : double);
var k : integer;
    sum, tmp : double;

begin
  if numberofhabitats > 0 then
  begin
    sum := 0;
    if dispersalpreference = 0 then sum := numberofhabitats - 1  // to simplify calculation if there is no dispersalpreference

    else for k := 0 to numberofhabitats-1 do // needs working on!!!!!!!!!!!!!!!!!!!
    begin
        tmp := Power((sqrdistance(i,j,habitats[k,0],habitats[k,1])), - 0.5 * dispersalpreference);
        if tmp > 0 then sum := sum + 1/tmp;
        //if sqrdistance > homerange then sum := sum + 1/ (abs(populatedhabitats[l,0]-i)-1, abs(populatedhabitats[l,1]-j)]
    end;

    sum := 1/sum * disprate;
    for k := 0 to numberofhabitats-1 do
    begin
        if (i <> habitats[k,0]) and (j <> habitats[k,1]) then
        ptemp[habitats[k,0],habitats[k,1]] := ptemp[habitats[k,0],habitats[k,1]]
        + sum * distancematrix[abs(unpopulatedhabitats[k,0]-i), abs(unpopulatedhabitats[k,1]-j)];
    end;
  end;
end;


procedure updatePopulation;
// one time step of the ecological model
var i,j,k: integer;
    presentmortality:double;
    lognormalmu, lognormalsigma: double; // input for the gaussian function to produce lognormal distr.
    sum : double;

begin

  // lognormal distribution for the mortality, calculates first the right parameters for the use of
  // the built in randg-function

  {lognormalmu:= ln(mortality)- 0.5* ln(1+rho/sqr(mortality));
  lognormalsigma := sqrt(ln(1+rho/sqr(mortality)));
  presentmortality := tanh(exp(randg(lognormalmu, lognormalsigma))) ;

  mortalitystatistics[ecocounter]:= presentmortality ;

  begin
  for i:=1 to dim do for j:=1 to dim do
      begin
      cneighbour[i,j]:=0;
      if i-1>=1 then ii:=i-1 else ii:=dim; cneighbour[i,j]:=cneighbour[i,j]+mortalitymatrix[ii,j]; cneighbour[ii,j]:=cneighbour[ii,j]+mortalitymatrix[i,j];
      if j-1>=1 then jj:=j-1 else jj:=dim; cneighbour[i,j]:=cneighbour[i,j]+mortalitymatrix[ii,jj]; cneighbour[ii,jj]:=cneighbour[ii,jj]+mortalitymatrix[i,j];
      if j+1<=dim then jj:=j+1 else jj:=1; cneighbour[i,j]:=cneighbour[i,j]+mortalitymatrix[ii,jj]; cneighbour[ii,jj]:=cneighbour[ii,jj]+mortalitymatrix[i,j];
      if j-1>=1 then jj:=j-1 else jj:=dim; cneighbour[i,j]:=cneighbour[i,j]+mortalitymatrix[i,jj];  cneighbour[i,jj]:=cneighbour[i,jj]+mortalitymatrix[i,j];
      end;
  for i:=1 to dim do for j:=1 to dim do mortalitymatrix[i,j]:= (2*random*presentmortality + 0.5*cneighbour[i,j]/8 + 0.5*mortalitymatrix[i,j]);
  end;}


  // destroyed habitats cause dispersal!
  // entweder hier oder nach dem Dispersal sterben lassen!

  createPopulationList;

  if numberofpopulatedhabitats > 0 then
  begin
    createHabitatList;

    // MORTALITY
    mortalitycircles(2);
    //picture_draw;

    {for l:= 0 to (numberofpopulatedhabitats-1) do
    begin
      i:= populatedhabitats[l,0];j:= populatedhabitats[l,1];
      case x[i,j] of
         0: p[i,j] := 0;
         1: if random < presentmortality then p[i,j] := 0 ;
         else if showwarnings then showmessage('Update Pop : Falscher Wert in p - matrix');
      end;
    end; }

    createPopulationList;

    // dispersal of the destroyed patches

    for i:=1 to dim do  for j:=1 to dim do
    begin
      ptemp[i,j]:= 0;
      if (x[i,j]= 0) and (p[i,j]= 1) then
      begin
        disperse(i,j,1);
        p[i,j]:= 0;
      end;
    end;

    createPopulationList;

    // Dispersal

    for k:= 0 to (numberofpopulatedhabitats-1) do
    begin
      i:= populatedhabitats[k,0];j:= populatedhabitats[k,1];  // read the coordinates of a habitat
      disperse(i,j, dispersalrate );
    end;
    for k:= 0 to (numberofunpopulatedhabitats-1) do
    begin
      i:= unpopulatedhabitats[k,0];j:= unpopulatedhabitats[k,1];
      if ptemp[i,j] > random then p[i,j] := 1 ;
    end;
    // Statistics

    createPopulationList;
  end;
  if showvalues then populationstatistics[ecocounter]:= numberofpopulatedhabitats/xtot;
  ecocounter := ecocounter +1 ;
end;


// GAME THEORY
// ----------------------------------------------


procedure constants_calculateGT;
var i,j : integer;
    halfdim: integer;
begin
  sqrdim := sqr(dim);
  {for i := 0 to 8 do
  begin
    benefitarray[i]:= ecological_benefit(i)
  end; }
  desiredBenefit := lambda*(1+ecological_benefit(8));
  halfdim := floor(dim/2+0.01);
  //setlength(distancematrix, halfdim, halfdim);
  for i := 0 to maxdim do for j := 0 to maxdim do distancematrix[i,j]:= 0;

  for i := 0 to (halfdim) do for j := 0 to (halfdim) do
  begin
    distancematrix[i,j]:= dispersalProbability(sqrt(i*i+j*j))
  end;

  for i := (dim) downto (halfdim+1) do for j := 0 to (halfdim) do
  begin
    distancematrix[i,j]:= distancematrix[(dim-i),j]
  end;
  for i := 0 to (halfdim) do for j := (dim) downto (halfdim+1) do
  begin
    distancematrix[i,j]:= distancematrix[i,(dim-j)]
  end;
  for i := (dim) downto (halfdim+1) do for j := (dim) downto (halfdim+1) do
  begin
    distancematrix[i,j]:= distancematrix[(dim-i),(dim-j)]
  end;
end;

procedure initializePopulationGT;
var i,j : integer;
    fractionofpopulatedhabitats: real;
    //minimumneighbours  : integer;
begin
  // variables
  fractionofpopulatedhabitats:=0.4;
  //minimumneighbours:= 3 ;
  for i:=1 to dim do  for j:=1 to dim do
  begin
  if (x[i,j]= 1) and ( random < fractionofpopulatedhabitats) {and (m[i,j]>= minimumneighbours)}
    then
    begin
      if random < 0.01 then p[i,j]:= 1 else p[i,j]:= 2;
    end
    else p[i,j]:= 0
  end;
end;

procedure createHabitatListGT;
var i,j,k,l : integer;
begin
  k:=0; l:=0;
  for i:=1 to dim do  for j:=1 to dim do
  begin
    if x[i,j]=1 then
    begin
      habitats[k,0]:= i;
      habitats[k,1]:= j;
      k:= k + 1;
    end;
    if (p[i,j]=1) or (p[i,j]=2)  then
    begin
      populatedhabitats[l,0]:=i;
      populatedhabitats[l,1]:=j;
      l:=l+1;
    end;
  end;
  if (k <> xtot) and showwarnings then showmessage('Overflow in createHabitatList');
  numberofpopulatedhabitats := l;
end;


procedure createPopulatedHabitatListGT;
var i,j,l : integer;
begin
  l:=0;
  for i:=1 to dim do  for j:=1 to dim do
  begin
    if (p[i,j]=1) or (p[i,j]=2)  then
    begin
      populatedhabitats[l,0]:=i;
      populatedhabitats[l,1]:=j;
      l:=l+1;
    end;
  end;
  numberofpopulatedhabitats := l;
end;


procedure updatePopulationGT;
var i,j,k,l,r: integer;
    presentmortality:double;
    lognormalmu, lognormalsigma: double; // input for the gaussian function to produce lognormal distr.
    sum : array [0..2] of double;
    mfitness,mstrategy, strategy : array [0..2] of double;
    fitness ,help: double;
    mutationrate : double;

begin

  mutationrate := 0.005;

  // lognormal distribution for the mortality, calculates first the right parameters for the use of
  // the built in randg-function


  //lognormalmu:= ln(mortality)- 0.5* ln(1+rho/sqr(mortality));
  //lognormalsigma := sqrt(ln(1+rho/sqr(mortality)));
 // presentmortality := tanh(exp(randg(lognormalmu, lognormalsigma))) ;

  presentmortality := mortality ;


  // MORTALITY

  createHabitatListGT;
  for r := 0 to 2 do
  begin
    mstrategy[r]:= 0;
    mfitness[r]:=0;
  end;

  // calculate fitness

  for k:= 0 to (numberofpopulatedhabitats-1) do
  begin
    sum[0] := 0; sum[2] :=0; sum[1] := 0;
    i:= populatedhabitats[k,0];j:= populatedhabitats[k,1];
    sum[ p[i,j]] := sum[ p[i,j]] - distancematrix[0,0];
    for l:= 0 to (xtot-1) do
    sum[ p[habitats[l,0],habitats[l,1]] ] :=
    sum[ p[habitats[l,0],habitats[l,1]] ] + distancematrix[abs(habitats[l,0]-i), abs(habitats[l,1]-j)] ;
    help := 1/(sum[1]+sum[2]+sum[0]) ;
    for r := 0 to 2 do
    begin
    strategy[r] := sum[r]*help;
    mstrategy[r]:= mstrategy[r]+ strategy[r];
    end;
    if p[i,j] = 1 then fitness:=  (cc * strategy[1] + cd * strategy[2]) else
    if p[i,j] = 2 then fitness:=  (dc * strategy[1] + dd * strategy[2]);
    fitness := fitness + dd * strategy[0] ;
    fitnessmatrix[i,j]:= fitness;
    mfitness[p[i,j]]:= mfitness[p[i,j]] + fitness;
  end;

  // calculate patches which die

  for k:= 0 to (numberofpopulatedhabitats-1) do
  begin
  if fitnessmatrix[populatedhabitats[k,0],populatedhabitats[k,1]] * random  < presentmortality then
  p[populatedhabitats[k,0],populatedhabitats[k,1]] := 0 ;
  end;
  for r := 0 to 2 do
  begin
    if numberofpopulatedhabitats <> 0 then mstrategy[r]:= mstrategy[r]/numberofpopulatedhabitats else mstrategy[r]:=0 ;
    if mstrategy[r] <> 0 then mfitness[r]:=mfitness[r]/numberofpopulatedhabitats/mstrategy[r] else mfitness[r]:= 0 ;
  end;


  // POPULATION

  createpopulatedHabitatListGT;

  for k:= 0 to (xtot-1) do
  begin
    i:= habitats[k,0];j:= habitats[k,1];  // read the coordinates of a habitat
    if p[i,j] = 0  then      // if it is not populated, calculate the colonization probability
    begin
      sum[0] := 0; sum[2] :=0; sum[1] := 0;
      for l:= 0 to (numberofpopulatedhabitats-1) do
      begin
      sum[ p[populatedhabitats[l,0],populatedhabitats[l,1] ] ] := sum[ p[populatedhabitats[l,0],populatedhabitats[l,1]] ] +
      distancematrix[ (abs(j-populatedhabitats[l,1])),(abs(i-populatedhabitats[l,0]))] * intpower(fitnessmatrix[populatedhabitats[l,0],populatedhabitats[l,1]],1);
      end;
      for r:= 0 to 2 do if sum[r] < 0 then sum[r] := 0;
      if (sum[1]+sum[2]) <> 0 then
      begin
      strategy[2] := sum[2]/(sum[1]+sum[2]);
      strategy[1] := sum[1]/(sum[1]+sum[2]) ;
      if sum[1]+sum[2] > random then if strategy[2] > random + (2-4*random)*mutationrate then p[i,j]:= 2 else p[i,j] := 1;
      end;
      // random colonization
      // if (p[i,j] = 0) and (random < 0.0001) then if mutationrate < random then p[i,j]:= 2 else p[i,j] := 1;
    end;
  end;
end;







// MAIN -------------------------------------------------------------

procedure run(cos:integer);
// implements one total time step incl ecological and economic model
var i,j,k : integer;
var mx : integer;
begin
  loopsum := 0 ; mprice := 0; mx := 0; turncount:=0; meantotalcosts :=0; mprecision :=0;
  basic_calculations_do;
  for k:=1 to maxt do
  begin
    // economic Model

    if economicmodelon then
    begin
      for i:=1 to dim do  for j:=1 to dim do xold[i,j]:=x[i,j];
      randomize_costs(cos);
      if optimizationon then optimaltrade else
      begin
        for j:=1 to trades do
        begin
          trade;
          loopsum := loopsum + loops;
          count:=count+1 ;
        end;
      end;
      for i:=1 to dim do  for j:=1 to dim do if x[i,j]<>xold[i,j] then
      begin
        turncount := turncount +1;
      end;
      calculations_do;
      mprice := mprice + price;
      mprecision := mprecision + precision;
      mx := mx + xtot;
      meantotalcosts := meantotalcosts + totalcosts;
    end;
    // ecological model

    if ecologicalmodelon then
    begin
      if not economicmodelon then basic_calculations_do;
      updatePopulation;
      // noch was machen, dass die # habitate richtig angezeigt wird
    end;
  end;
  // statistics
  if economicmodelon then
  begin
    mprice := mprice/maxt;
    mxtot := mx/maxt;
    mprecision := mprecision/maxt;
    meantotalcosts := meantotalcosts / maxt;
    turnoverrate:= turncount/maxt/mxtot/2;
  end
  else
  begin
    mxtot := 1;
    for i:=1 to dim do  for j:=1 to dim do
    mxtot := mxtot + x[i,j];
  end;
end;

// -----------------------------------------------------------------
// EXPERIMENTS

procedure ecomodel;
// makes a PVA for the current setting
// erster Parameter : Anfangsverteilung
// zweiter Parameter : Kostenfunktion
var i,j,k,outputvariables: integer;
statistics : TextFile;
time : string;
begin
  time := FormatDateTime('yy-mm-dd_hh-mm-ss_',Now);
  ecologicalmodelon := true;
  write_ini('statistics\'+ time + 'popstat.ini');
  AssignFile(statistics, 'statistics\' + time + 'popstat.dat');
  rewrite(statistics);
  append(statistics);
  writeln(statistics, '# Remarks: population statistics' );
  writeln(statistics, '# samplesize', varytries );
// varies the value and runs the trade ---------------------------------------------------
  for j := 1 to varytries do
  begin
    constants_calculate;
    calculations_do;
    initializePopulation;
    initializeEcoStatistics(maxt*runs);
    ecocounter:=0;
    createHabitatList;
    for k := 1 to (runs) do run(costvariation);
    i:= 0;
    while populationstatistics[i]>0 do i := i + 1;
    writeln(statistics, j , chr(9), i);
  end;
  closefile(statistics) ;
end;





procedure einsdscan(rand, cos:integer);
// erster Parameter : Anfangsverteilung
// zweiter Parameter : Kostenfunktion
var i,j,k,l: integer;
    vst : real;
    point : ^double;

begin
  point :=@sigma ;
  //varyaround := point^ ;  // war unsinnig und sollte eigentlich nicht mehr gebraucht werden!
  point^ := varyaround - varyrange/2  ;
  vst := varyrange/varysteps   ;
  if showvalues then form1.runstatistics.Items.Clear;
  outfile_initialize;
  // varies the value and runs the trade ---------------------------------------------------
  for i := 1 to (varysteps +1)  do
  begin
    constants_calculate;
    if showvalues then parameters_show;
    for l:=1 to 5 do summary[l,1]:=0; for l:=1 to 5 do summary[l,2]:=0;
    for j := 1 to varytries do
    begin
      randomize_patches(rand);
      randomize_costs(cos);

      for k := 1 to (runs) do run(cos);
      run(cos);

      calculations_do;
      summary[1,1]:= summary[1,1] + turnoverrate   ;
      summary[2,1]:= summary[2,1] + totalcosts ;
      summary[3,1]:= summary[3,1] + connectivity  ;
      summary[4,1]:= summary[4,1] + mxtot/sqrdim  ;
      summary[5,1]:= summary[5,1] + mprice   ;
      summary[1,2]:= summary[1,2] + sqr(turnoverrate)  ;
      summary[2,2]:= summary[2,2] + sqr(totalcosts)     ;
      summary[3,2]:= summary[3,2] + sqr(connectivity ) ;
      summary[4,2]:= summary[4,2] + sqr(mxtot/sqrdim ) ;
      summary[5,2]:= summary[5,2] + sqr(mprice )  ;
    end;
    for l := 1 to 5 do summary[l,2]:= sqrt(abs(summary[l,2]-sqr(summary[l,1])/varytries) /(varytries-1));
    for l := 1 to 5 do summary[l,1]:= summary[l,1] /varytries    ;
    if showvalues then
    begin
      picture_draw(0,0);
      logwindow_show(i);
    end;
    outfile_write;
    point^ := point^ + vst ;
  end;
  closefile(outfile) ;
  point^:= varyaround;
  if showvalues then
  begin
    parameters_show;
    form1.RunStatistics.Items.insert(0,'Finished');
  end;
end;


procedure einsdscanoptimization(rand, cos:integer);
// erster Parameter : Anfangsverteilung
// zweiter Parameter : Kostenfunktion
var i,j,l, outputvariables: integer;
    vst : real;
    point : ^double;

begin
  point :=@sigma ;
  //varyaround := point^ ;  // war unsinnig und sollte eigentlich nicht mehr gebraucht werden!
  point^ := varyaround - varyrange  ;
  vst := varyrange/(varysteps -1)  ;
  outfile_initialize;
  outputvariables :=5;

  // varies the value and runs the trade ---------------------------------------------------
  for i := 1 to (varysteps )  do
  begin
    constants_calculate;
    if showvalues then parameters_show;
    for l:=1 to outputvariables do summary[l,1]:=0; for l:=1 to outputvariables do summary[l,2]:=0;
    for j := 1 to varytries do
    begin
      randomize_patches(rand);
      randomize_costs(cos);

      run(cos);

      calculations_do;
      summary[1,1]:= summary[1,1] + turnoverrate   ;
      summary[2,1]:= summary[2,1] + totalcosts ;
      summary[3,1]:= summary[3,1] + connectivity  ;
      summary[4,1]:= summary[4,1] + mxtot/sqrdim  ;

      summary[1,2]:= summary[1,2] + sqr(turnoverrate)  ;
      summary[2,2]:= summary[2,2] + sqr(totalcosts)     ;
      summary[3,2]:= summary[3,2] + sqr(connectivity ) ;
      summary[4,2]:= summary[4,2] + sqr(mxtot/sqrdim ) ;

    end;

    for l := 1 to outputvariables do
        begin
          summary[l,2]:= sqrt(abs(summary[l,2]-sqr(summary[l,1])/varytries) /(varytries-1));
          summary[l,1]:= summary[l,1] /varytries ;
    end;

    outfile_write;
    point^ := point^ + vst ;
  end;
  closefile(outfile) ;
  point^:= varyaround;
end;



procedure einsdscaneco(rand, cos:integer);
// erster Parameter : Anfangsverteilung
// zweiter Parameter : Kostenfunktion
var i,j,k,l,n,outputvariables: integer;
    ecorun : boolean;
    vst : real;
    point : ^double;

begin
  point :=@w ;
  //varyaround := point^ ;   // war unsinnig und sollte eigentlich nicht mehr gebraucht werden!
  point^ := varyaround - varyrange  ;
  vst := varyrange/(varysteps-1)*2   ;
  outfile_initialize;
  outputvariables := 7;
  ecorun :=  ecologicalmodelon;
  // varies the value and runs the trade ---------------------------------------------------
  for i := 1 to (varysteps )  do
  begin
    constants_calculate;
    if showvalues then parameters_show;
    for l:=1 to outputvariables do summary[l,1]:=0; for l:=1 to outputvariables do summary[l,2]:=0;
    for j := 1 to varytries do
    begin
      randomize_patches(startdistribution);
      randomize_costs(0);
      ecocounter:=0;

      ecologicalmodelon := false;
      for k := 1 to (runs) do run(costvariation);

      if ecorun then
      begin
        ecologicalmodelon := true;
        initializePopulation;
        ecocounter:=0;
        initializeEcoStatistics(maxt);
      end;

      run(costvariation);

      calculations_do;

      // assign values for mean output
      summary[1,1]:= summary[1,1] + turnoverrate   ;
      summary[2,1]:= summary[2,1] + totalcosts ;
      summary[3,1]:= summary[3,1] + connectivity  ;
      summary[4,1]:= summary[4,1] + mxtot/sqrdim  ;
      summary[5,1]:= summary[5,1] + mprice   ;
      summary[6,1]:= summary[6,1] + numberofpopulatedhabitats ;
      n:=0  ; repeat n := n+1 until populationstatistics[n]<0.000001; summary[7,1]:= summary[7,1] + n ;

      // sum squares for variance output

      summary[1,2]:= summary[1,2] + sqr(turnoverrate)   ;
      summary[2,2]:= summary[2,2] + sqr(totalcosts) ;
      summary[3,2]:= summary[3,2] + sqr(connectivity)  ;
      summary[4,2]:= summary[4,2] + sqr(mxtot/sqrdim)  ;
      summary[5,2]:= summary[5,2] + sqr(mprice)   ;
      summary[6,2]:= summary[6,2] + sqr(numberofpopulatedhabitats) ;
      summary[7,2]:= summary[7,2] + sqr(n) ;
    end;
    for l := 1 to outputvariables do summary[l,2]:= sqrt(abs(summary[l,2]-sqr(summary[l,1])/varytries) /(varytries-1));
    for l := 1 to outputvariables do summary[l,1]:= summary[l,1] /varytries    ;
    outfile_write;
    point^ := point^ + vst ;
  end;
  closefile(outfile) ;
  point^:= varyaround; ecologicalmodelon := ecorun;
end;




procedure einsdscanintensity(rand, cos:integer);   // scant verschieden Intesities
// erster Parameter : Anfangsverteilung
// zweiter Parameter : Kostenfunktion
var i,j,k,l: integer;
    vst : integer;
    point : ^integer;

begin
  point :=@intensity ;
  //varyaround := point^ ;   // war unsinnig und sollte eigentlich nicht mehr gebraucht werden!
  point^ := 1;
  vst := 1   ;
  if showvalues then form1.runstatistics.Items.Clear;
  outfile_initialize;
  // varies the value and runs the trade ---------------------------------------------------
  for i := 1 to (varysteps)  do
  begin
    constants_calculate;
    if showvalues then parameters_show;
    for l:=1 to 5 do summary[l,1]:=0; for l:=1 to 5 do summary[l,2]:=0;
    for j := 1 to varytries do
    begin
      randomize_patches(rand);
      randomize_costs(cos);

      for k := 1 to (runs) do run(cos);
      run(cos);

      calculations_do;
      summary[1,1]:= summary[1,1] + turnoverrate   ;
      summary[2,1]:= summary[2,1] + totalcosts ;
      summary[3,1]:= summary[3,1] + connectivity  ;
      summary[4,1]:= summary[4,1] + mxtot/sqrdim  ;
      summary[5,1]:= summary[5,1] + mprice   ;

      summary[1,2]:= summary[1,2] + sqr(turnoverrate)  ;
      summary[2,2]:= summary[2,2] + sqr(totalcosts)     ;
      summary[3,2]:= summary[3,2] + sqr(connectivity ) ;
      summary[4,2]:= summary[4,2] + sqr(mxtot/sqrdim ) ;
      summary[5,2]:= summary[5,2] + sqr(mprice )  ;
    end;
    for l := 1 to 5 do summary[l,2]:= sqrt(abs(summary[l,2]-sqr(summary[l,1])/varytries) /(varytries-1));
    for l := 1 to 5 do summary[l,1]:= summary[l,1] /varytries    ;
    if showvalues then
    begin
      picture_draw(0,0);
      logwindow_show(i);
    end;
    outfile_write;
    point^ := point^ + vst ;
  end;
  closefile(outfile) ;
  point^:= 1;
  if showvalues then
  begin
    parameters_show;
    form1.RunStatistics.Items.insert(0,'Finished');
  end;
end;


procedure zweidscan(rand, cos:integer);
// erster Parameter : Anfangsverteilung
// zweiter Parameter : Kostenfunktion
var i,j,k,l,m, outputvariables: integer;
    vst1,vst2 : real;
    point1 : ^double;
    point2 : ^double;

begin
  point1 :=@sigma ;
  point2 :=@w ;
  //varyaround2d1 := point1^ ;      // war unsinnig und sollte eigentlich nicht mehr gebraucht werden!
  //varyaround2d2 := point2^ ;      // war unsinnig und sollte eigentlich nicht mehr gebraucht werden!
  point1^ := varyaround2d1 - varyrange2d1  ;
  point2^ := varyaround2d2 - varyrange2d2  ;
  vst1 := varyrange2d1/(varysteps2d1-1)*2   ;
  vst2 := varyrange2d2/(varysteps2d2-1)*2   ;
  outfile_initialize;
  outputvariables := 5;
// varies the value and runs the trade ---------------------------------------------------
  for i := 1 to (varysteps2d1 )  do
  begin
    for m := 1 to (varysteps2d2 )  do
    begin
      constants_calculate;
      for l:=1 to outputvariables do summary[l,1]:=0; for l:=1 to outputvariables do summary[l,2]:=0;
      for j := 1 to varytries do
      begin
        randomize_patches(rand);
        randomize_costs(cos);

        basic_calculations_do;

        if not optimizationon then for k := 1 to (runs) do run(cos);
        run(cos);
        calculations_do;
        summary[1,1]:= summary[1,1] + turnoverrate   ;
        summary[2,1]:= summary[2,1] + totalcosts      ;
        summary[3,1]:= summary[3,1] + connectivity  ;
        summary[4,1]:= summary[4,1] + mxtot/sqrdim  ;
        summary[5,1]:= summary[5,1] + mprice   ;
        if not optimizationon then
        begin
        summary[1,2]:= summary[1,2] + sqr(turnoverrate)  ;
        summary[2,2]:= summary[2,2] + sqr(totalcosts)     ;
        summary[3,2]:= summary[3,2] + sqr(connectivity ) ;
        summary[4,2]:= summary[4,2] + sqr(mxtot/sqrdim ) ;
        summary[5,2]:= summary[5,2] + sqr(mprice )  ;
        end;
      end;
      if not optimizationon then for l := 1 to outputvariables do summary[l,2]:= sqrt(abs(summary[l,2]-sqr(summary[l,1])/varytries) /(varytries-1));
      for l := 1 to outputvariables do summary[l,1]:= summary[l,1] /varytries ;

      outfile_write;
      point1^ := point1^ + vst1 ;
    end;
    point1^ := varyaround2d1 - varyrange2d1  ;
    point2^ := point2^ + vst2 ;
  end;
  closefile(outfile) ;
  point1^:= varyaround2d1;
  point2^:= varyaround2d2;
end;


procedure optincentivescan(rand, cos:integer);
// erster Parameter : Anfangsverteilung
// zweiter Parameter : Kostenfunktion
var i,j,k,l,m,n,nn, red, outputvariables,o: integer;
    vst1,vst2, vst3 : real;
    point1 : ^double;
    point2 : ^double;
    point3 : ^double;
begin
  point1 :=@lambda;
  point2 :=@w ;
  point3 :=@interactionlength ;
  if (varytries mod reduce) <> 0 then application.terminate; 
  if varysteps > 1 then point1^ := varyaround - varyrange else point1^ := varyaround ;
  if varysteps2d1 > 1 then point2^ := varyaround2d1 - varyrange2d1 else point2^ := varyaround2d1 ;
  if varysteps2d2 > 1 then point3^ := varyaround2d2 - varyrange2d2 else point3^ := varyaround2d2 ;
  if varysteps > 1 then vst1 := 2*varyrange/(varysteps-1) else vst1 := 0 ;
  if varysteps2d1 > 1 then vst2 := 2*varyrange2d1/(varysteps2d1-1)  else vst2 := 0  ;
  if varysteps2d2 > 1 then vst3 := 2*varyrange2d2/(varysteps2d2-1)  else vst3 := 0  ;
  outputvariables := 7;
  outfile_initialize;
// varies the value and runs the trade ---------------------------------------------------
  for n := 1 to (varysteps )  do
  begin
    for i := 1 to (varysteps2d1 )  do
    begin
      for m := 1 to (varysteps2d2 )  do
      begin

        constants_calculate;
        for l:=1 to outputvariables do summary[l,1]:=0; for l:=1 to outputvariables do summary[l,2]:=0;
        // get statistics on the parameter combinations

        for j := 1 to round(varytries / reduce) do
        begin
          reset;
          randomize_patches(rand);
          randomize_costs(cos);
          basic_calculations_do;
          ecologicalmodelon := false;
          for k := 1 to (runs) do run(cos);
          ecologicalmodelon := true;
          for red := 1 to reduce do
          begin
            initializePopulation;
            initializeEcoStatistics(maxt);
            ecocounter:=0;
            createHabitatList;

            run(cos);
            calculations_do;
            summary[1,1]:= summary[1,1] + turnoverrate   ;
            summary[2,1]:= summary[2,1] + totalcosts      ;
            summary[3,1]:= summary[3,1] + connectivity  ;
            summary[4,1]:= summary[4,1] + mxtot/sqrdim  ;
            summary[5,1]:= summary[5,1] + mprice   ;
            summary[6,1]:= summary[6,1] + numberofpopulatedhabitats/(mxtot+0.0000000001) ;
            if numberofpopulatedhabitats > 0 then summary[7,1]:= summary[7,1]+1;

            summary[1,2]:= summary[1,2] + sqr(turnoverrate)  ;
            summary[2,2]:= summary[2,2] + sqr(totalcosts)     ;
            summary[3,2]:= summary[3,2] + sqr(connectivity ) ;
            summary[4,2]:= summary[4,2] + sqr(mxtot/sqrdim ) ;
            summary[5,2]:= summary[5,2] + sqr(mprice )  ;
            summary[6,2]:= summary[6,2] + sqr(numberofpopulatedhabitats/(mxtot+0.0000000001));
            summary[7,2]:= 0
          end;  
        end;

        for l := 1 to outputvariables do
        begin
          summary[l,2]:= sqrt(abs(summary[l,2]-sqr(summary[l,1])/varytries) /(varytries-1));
          summary[l,1]:= summary[l,1] /varytries ;
        end;
        outfile_write;
        point3^ := point3^ + vst3 ;
      end;
      if varysteps2d2 > 1 then point3^ := varyaround2d2 - varyrange2d2 ;
      point2^ := point2^ + vst2 ;
    end;
    if varysteps2d2 > 1 then point3^ := varyaround2d2 - varyrange2d2 ;
    if varysteps2d1 > 1 then point2^ := varyaround2d1 - varyrange2d1  ;
    point1^ := point1^ + vst1 ;
  end;
  closefile(outfile) ;
  point2^:= varyaround2d1;
  point3^:= varyaround2d2;
  point1^:= varyaround ;
end;



//Buttons   -------------------------------------------------



procedure TForm1.Button1Click(Sender: TObject);  // button distribution reset
begin
  reset;
end;


procedure TForm1.ButtonEcologicalModelResetClick(Sender: TObject);    // button eco distribution reset
begin
  if ecologicalmodelon then
  begin
     constants_calculate;
     calculations_do;
     initializePopulation;
     createHabitatList;
     parameters_show;
     values_show;
     picture_draw(0,0);
  end
  else showmessage('Ecological Model is not on')
end;





procedure TForm1.Button2Click(Sender: TObject);   // button run
var k,i,j,l: integer;
    arr : array [0..300] of integer;
begin
  breakrun := false;
  logwindow_initialize;
  constants_calculate;
  parameters_show;
  calculations_do;
  form1.RunStatistics.Items.insert(0,'Starttime: '+TimeTOStr(Now));
  storetime := Now;
  basic_calculations_do;
  count:= 0;
  if CheckBoxCalculateCostDistribution.checked then form1.series3.clear;
  if showvalues then
  begin
  ecocounter:=0;  initializeEcoStatistics(maxt*runs);
  end;
  k := 0;
  while (k < runs) and not breakrun do
  begin
    if k = runs then calculations := true else calculations := false;
    run(form1.RadioGroupCosts.Itemindex);
    picture_draw(k,0);
    calculations_do;
    values_show;
    logwindow_show(k);
    k := k+1;
    //Application.ProcessMessages;
  end;
  if showvalues then
  begin
    form1.RunStatistics.Items.insert(0,'Finished after: ' + FormatDateTime('hh:nn:ss:zzz', ( Now-storetime)) + 'ms');


    if CheckBoxCalculateCostDistribution.checked then
    begin
      form1.series3.clear;
      for i := 0 to 300 do  arr[i]:=0;
      for i := 1 to dim do for j := 1 to dim do
      begin
         l := round((1-c[i,j])*150/strtofloat(form1.editcostdistributionrange.text));
         if l < -150 then l := -150;
         if l > 150 then l:=150 ;
         arr[l+150]:= arr[l+150] + 1
      end;
      for l := 0 to 300 do with form1.Series3 do AddXY((l-150)/150*strtofloat(form1.editcostdistributionrange.text),arr[l]);
    end;

    if ecologicalmodelon then
    begin
      form1.series1.clear;
      form1.series2.clear;
      for i:=0 to maxt*runs -1  do
      begin
        if form1.checkboxshowmortality.checked then with form1.Series1 do
            AddXY(i,mortalitystatistics[i]);
        if form1.checkboxshowpopulation.checked then with form1.Series2 do
            AddXY(i,populationstatistics[i]);
      end;
    end
  end;
end;





procedure TForm1.varyRunClick(Sender: TObject);  // button vary Run for variation of parameters
begin
   case variationchoice of
   0: einsdscan(StartDistribution,costvariation);
   1: zweidscan(StartDistribution,costvariation);
   2: optincentivescan(StartDistribution,costvariation);
   3: einsdscaneco(StartDistribution,costvariation);
   end;
end;






// ------------------
// INI Dialogs
// --------------------

procedure TForm1.form1close(Sender: TObject; var Action: TCloseAction);
begin
inifilename := ExtractFilePath(ParamStr(0))+ 'Project1.ini';
write_ini(inifilename);
end;

procedure TForm1.ButtonOpenIniClick(Sender: TObject);
var
  openDialog : TOpenDialog;    // Open dialog variable
begin
  try
    OpenDialog:= TOpenDialog.Create(self);
    OpenDialog.InitialDir := ExtractFileDir(inifilename);
    openDialog.Filter :='Ini Files|*.ini|All files|*.*';
    // Display the open file dialog
    if OpenDialog.Execute then
    begin
      inifilename := OpenDialog.FileName;
    end;
  finally
  OpenDialog.Free;
  end;
  read_ini;
  parameters_show;
  values_show;
end;

procedure TForm1.Button3Click(Sender: TObject);  // read ini
begin
  read_ini;
  parameters_show;
  values_show;
end;


{procedureTform1.Speichern1.CLick(Sender:TObject);
begin
if SaveDialog1.Execute then begin
 Memo1.Lines.SaveToFile(SaveDialog1.FileName);
 end;
end;
 }


procedure TForm1.ButtonSaveIniClick(Sender: TObject);   // write ini
begin
write_ini(inifilename);
parameters_show;
values_show;
end;

procedure TForm1.ButtonIniSaveAsClick(Sender: TObject);
var
  saveDialog : TSaveDialog;    // save dialog variable
begin
  try
    saveDialog:= TSaveDialog.Create(self);
    saveDialog.InitialDir := ExtractFileDir(inifilename);
    saveDialog.Filter :='Ini Files|*.ini|All files|*.*';
    saveDialog.DefaultExt := '.ini';
    // Display the save file dialog
    if saveDialog.Execute then
    begin
      inifilename := saveDialog.FileName;
    end;
  finally
  saveDialog.Free;
  end;
  write_ini(inifilename);
  parameters_show;
  values_show;
end;


procedure TForm1.ButtonBreakClick(Sender: TObject);
begin
breakrun := true;
end;

procedure TForm1.ButtonLogWindowClearClick(Sender: TObject);
begin
 form1.RunStatistics.Clear;
end;



//----------------------------------------------------------------------------------
//  Edits
//---------------------------------------------------------------------------------

procedure eingabe_error ;
begin
   if showwarnings then form1.RunStatistics.Items.insert(0, 'String error in the changed field');
end;

procedure warning(code:integer);
begin
  if showwarnings then
  begin
    case code of
      1: form1.RunStatistics.Items.insert(0, 'String error in the changed field');

      else showmessage('Funktion showwarnings aufgerufen mit falschem Fehlercode')
    end;
  end;
end;



procedure enableForms;
begin
    if showvalues then
    begin
      if (serial = 0) or (serial = 2) then form1.editintensity.enabled:= false
      else form1.editintensity.enabled:= true ;
      if altruismon = true then form1.editaltruism.enabled := true else form1.editaltruism.enabled := false;

      if (costvariation = 4) then form1.editcostsize.enabled := true else form1.editcostsize.enabled := false;
      if (costvariation = 2) or (costvariation = 1) then form1.edit11.Enabled := true else form1.edit11.enabled := false;
      case tradeconstraints of
      0:
        begin
          form1.label9.caption:=' % of Habitats if clumped';
          form1.label4.caption:='Desired Ecopoints';
        end;
      1:
        begin;
          form1.label9.caption:='Desired Totalcosts';
          form1.label4.caption:='Desired Totalcosts';
        end;
      end;
      case costvariation of
      1:
        begin;
          form1.label2.caption:='Max Range';
          form1.label11.caption:='Step length';
        end;
      2:
        begin;
          form1.label2.caption:='Step length';
          form1.label11.caption:='Cost Correlation Strength';
        end;
      else
        begin
          form1.label2.caption:='Cost Variation Range';
          form1.label11.caption:='Disabled';
        end;

      end;


      form1.editmortality.enabled:= form1.checkboxecologicalmodelon.checked;
      form1.editmortalitycorrelation.enabled:=  form1.checkboxecologicalmodelon.checked;
      form1.editdispersaldistance.enabled:=  form1.checkboxecologicalmodelon.checked;
      form1.editdispersalrate.enabled:=  form1.checkboxecologicalmodelon.checked;
      case variationchoice of
      0: form1.LabelVariationParameter.Caption:= 'Parameter1                inactive               inactive';
      1: form1.LabelVariationParameter.Caption:= 'inactive                Parameter1               Parameter2';
      2: form1.LabelVariationParameter.Caption:= 'goal labda               strength w               length';
      3: form1.LabelVariationParameter.Caption:= 'Parameter1                inactive               inactive';
      end;
   end;
end;

procedure TForm1.varound2d1Change(Sender: TObject);
var code: integer;
begin
  val(varound2d1.text,varyaround2d1,code);
  if code <> 0 then eingabe_error ;
end;

procedure TForm1.varound2d2Change(Sender: TObject);
var code: integer;
begin
  val(varound2d2.text,varyaround2d2,code);
  if code <> 0 then eingabe_error ;
end;


procedure TForm1.CheckBoxOptimizationShowClick(Sender: TObject);
begin
   optimizationshow := Form1.CheckBoxOptimizationShow.checked;
end;

procedure TForm1.EditStartTemperatureChange(Sender: TObject);
var code: integer;
begin
  val(EditStartTemperature.text,StartTemperature,code);
  if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditTemperatureDecayChange(Sender: TObject);
var code: integer;
begin
  val(EditTemperatureDecay.text,TemperatureDecay,code);
  if code <> 0 then eingabe_error ;
end;



procedure TForm1.Chart1Click(Sender: TObject);
begin
  chart1.copytoclipboardbitmap ;
end;

procedure TForm1.Image1Click(Sender: TObject);
begin
  clipboard.assign(form1.image1.picture);
end;

procedure TForm1.EditDispersalPreferenceChange(Sender: TObject);
var code: integer;
begin
  val(Editdispersalpreference.text,dispersalpreference,code);
  if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditvradiusChange(Sender: TObject);
var code: integer;
begin
  val(Editvradius.text,vradius,code);
  if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditGTccChange(Sender: TObject);
var code: integer;
begin
  val(EditGTcc.text,cc,code);
  if code <> 0 then eingabe_error ;
end;


procedure TForm1.EditGTcdChange(Sender: TObject);
var code: integer;
begin
 val(EditGTcd.text,cd,code);
 if code <> 0 then eingabe_error ;
end;


procedure TForm1.EditGTdcChange(Sender: TObject);
var code: integer;
begin
  val(EditGTdc.text,dc,code);
  if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditGTddChange(Sender: TObject);
var code: integer;
begin
  val(EditGTdd.text,dd,code);
  if code <> 0 then eingabe_error ;
end;


procedure TForm1.Edit7Change(Sender: TObject);     // input dimension
var code: integer;
begin
  val(edit7.text,dim,code);
  if code <> 0 then eingabe_error ;
end;

procedure TForm1.Edit1Change(Sender: TObject); // input w
var code: integer;
begin
   val(edit1.text,w,code);
   if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditLengthChange(Sender: TObject); // input length
var code: integer;
begin
   val(editLength.text,interactionlength,code);
   if code <> 0 then eingabe_error;
   constants_calculate;
end;

procedure TForm1.Edit2Change(Sender: TObject);  // input sigma
var code: integer;
begin
   val(edit2.text,sigma,code);
   if code <> 0 then eingabe_error ;
end;

procedure TForm1.Edit9Change(Sender: TObject);    // input Lambda
var code: integer;
begin
   val(edit9.text,lambda,code);
   if code <> 0 then eingabe_error ;
end;

procedure TForm1.Edit11Change(Sender: TObject);    // input pricestability
var code: integer;
begin
  val(edit11.text,pricestability,code);
  if code <> 0 then eingabe_error ;
end;

procedure TForm1.TradeRunsChange(Sender: TObject);     // input trades
var code: integer;
begin
    val(traderuns.text,trades,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.Edit3Change(Sender: TObject);   // input maxt (intermediate runs)
var code: integer;
begin
   val(edit3.text,maxt,code);
   if code <> 0 then eingabe_error ;
end;

procedure TForm1.Edit10Change(Sender: TObject);  //input runs
var code: integer;
begin
   val(edit10.text,runs,code);
   if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditIntensityChange(Sender: TObject);      //intesity
var code: integer;
begin
    val(editintensity.text,intensity,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditPrecisionChange(Sender: TObject);
var code: integer;
begin
    val(editprecision.text,precision,code);
    if code <> 0 then eingabe_error ;
end;



//graphics


procedure TForm1.EditZoomChange(Sender: TObject);
var code: integer;
begin
    val(editzoom.text,zoom,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.UpDownZoomClick(Sender: TObject; Button: TUDBtnType);
begin

end;


//options

procedure TForm1.CheckBoxShowWarningsClick(Sender: TObject);
begin
if Form1.CheckBoxShowWarnings.Checked then showwarnings := true
else showwarnings := false;
end;

procedure TForm1.RadioGroupSerialClick(Sender: TObject);
begin
  serial:= form1.RadioGroupSerial.Itemindex;
  enableforms;
end;

procedure TForm1.RadioGroupVariationChoiceClick(Sender: TObject);
begin
  variationchoice := RadioGroupVariationChoice.itemindex;
  enableforms;
end;



procedure TForm1.CheckBoxEcologicalModelOnClick(Sender: TObject);
begin
  if Form1.CheckBoxEcologicalmodelon.Checked then ecologicalmodelon := true
  else ecologicalmodelon := false;
  enableforms;
end;

procedure TForm1.CheckBoxEconomicModelOnClick(Sender: TObject);
begin
  if Form1.CheckBoxEconomicmodelon.Checked then economicmodelon := true
  else economicmodelon := false;
end;

procedure TForm1.RadioGroupCostsClick(Sender: TObject);
begin
costvariation:= Form1.RadioGroupCosts.itemindex;
enableforms;
end;

procedure TForm1.RadioGroupStartDistributionClick(Sender: TObject);
begin
StartDistribution:= Form1.RadioGroupStartDistribution.itemindex;
end;

procedure TForm1.Memo1Change(Sender: TObject);
begin
 remarks := Form1.Memo1.Text;
end;



//ecological model

procedure TForm1.EditMortalityChange(Sender: TObject);
var code: integer;
begin
    val(editmortality.text,mortality,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditMortalityCorrelationChange(Sender: TObject);
var code: integer;
begin
    val(editmortalitycorrelation.text,rho,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditDispersalRateChange(Sender: TObject);
var code: integer;
begin
    val(editdispersalrate.text,dispersalrate,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditDispersalDistanceChange(Sender: TObject);
var code: integer;
begin
    val(editdispersaldistance.text,alpha,code);
    if code <> 0 then eingabe_error ;
end;

//---------variation--------------------

procedure TForm1.vRangeChange(Sender: TObject);   // input variation range
var code: integer;
begin
    val(vrange.text,varyrange,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.vStepsChange(Sender: TObject);   // input variation steps
var code: integer;
begin
    val(vsteps.text,varysteps,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.variationTriesChange(Sender: TObject);  //input variation tries
var code: integer;
begin
    val(variationTries.text,varytries,code);
    if varytries < 2 then
    begin
       varytries := 2;
       form1.RunStatistics.Items.insert(0,'varitries adjusted, must not be smaller 2');
       parameters_show;
    end;
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.vAroundChange(Sender: TObject);
var code: integer;
begin
    val(vAround.text,varyaround,code);
    if code <> 0 then eingabe_error ;
end;

//--  variation 2d -------------------------

procedure TForm1.vRange2d1Change(Sender: TObject);   // input variation range
var code: integer;
begin
    val(vrange2d1.text,varyrange2d1,code);
    if code <> 0 then eingabe_error ;

end;

procedure TForm1.vSteps2d1Change(Sender: TObject);   // input variation steps
var code: integer;
begin
    val(vsteps2d1.text,varysteps2d1,code);
    if code <> 0 then  eingabe_error ;
end;


procedure TForm1.vRange2d2Change(Sender: TObject);   // input variation range
var code: integer;
begin
    val(vrange2d2.text,varyrange2d2,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.vSteps2d2Change(Sender: TObject);   // input variation steps
var code: integer;
begin
    val(vsteps2d2.text,varysteps2d2,code);
    if code <> 0 then eingabe_error ;
end;



procedure TForm1.FilenameChange(Sender: TObject);    // input filename
begin
statisticsfilename := form1.filename.text
end;


procedure TForm1.RadioGroupShowValuesClick(Sender: TObject);
begin
if form1.RadioGroupShowValues.Itemindex = 0 then showvalues := true
else showvalues := false;
end;

procedure TForm1.EditAltruismChange(Sender: TObject);
var code: integer;
begin
    val(editaltruism.text,altruism,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.CheckBoxAltruismonClick(Sender: TObject);
begin
  altruismon := Form1.CheckBoxAltruismon.checked;
  enableforms;
end;

procedure TForm1.CheckBoxOptimizationOnClick(Sender: TObject);
begin
   optimizationon := Form1.CheckBoxoptimizationon.checked;
end;


procedure TForm1.EditOptimizationPercentageChange(Sender: TObject);
var code: integer;
begin
    val(EditOptimizationPercentage.text,OptimizationPercentage,code);
    if code <> 0 then eingabe_error ;
end;


procedure TForm1.RadioGrouptradeconstraintsClick(Sender: TObject);
begin
  tradeconstraints := RadioGrouptradeconstraints.Itemindex;
  enableforms;
end;


procedure TForm1.CheckBoxIndependentrunonClick(Sender: TObject);
begin
  independentrunon := CheckBoxIndependentrunon.checked
end;




procedure TForm1.n1stvariationClick(Sender: TObject);
begin
  variation1:= Tmenuitem(Sender).Tag;
  enableforms;
end;

procedure TForm1.Popupmenuvariation2itemClick(Sender: TObject);
begin
  variation2:= Tmenuitem(Sender).Tag;
  enableforms;
end;

procedure TForm1.EditCostSizeChange(Sender: TObject);
var code: integer;
begin
    val(EditCostSize.text,Costsize,code);
    if code <> 0 then eingabe_error ;
end;

procedure TForm1.EditmradiusChange(Sender: TObject);
var code: integer;
begin
  val(Editmradius.text,mradius,code);
  if code <> 0 then eingabe_error ;
end;


procedure TForm1.EditReduceChange(Sender: TObject);
var code: integer;
begin
  val(EditReduce.text,reduce,code);
  if code <> 0 then eingabe_error ;
end;


// Game Theory -------------------------------------------------------

procedure TForm1.ButtonGTResetClick(Sender: TObject);
begin
   initializePopulationGT;
   parameters_show;
   picture_draw(0,0);
end;

procedure TForm1.ButtonGTRunClick(Sender: TObject);
var k,i: integer;
begin
  basic_calculations_do;
  constants_calculateGT;
  parameters_show;
  count:= 0;
  for k := 1 to runs do
  begin
    loopsum := 0 ; mprice := 0; turncount:=0; meantotalcosts :=0;
    for i:=1 to maxt do updatePopulationGT;
    picture_draw(0,0);
    //values_show;
  end;
end;

procedure TForm1.Button4Click(Sender: TObject);
begin
  costfield_write;
end;

procedure TForm1.Button5Click(Sender: TObject);
begin
  ecomodel;
end;


procedure TForm1.Button6Click(Sender: TObject);
begin
   costseries_write
end;





//main-------------------------------------------------------------

procedure reset;
begin
   constants_calculate;
   randomize_costs(0);
   randomize_patches(startdistribution) ;
   calculations_do;
   parameters_show;
   values_show;
   picture_draw(0,0);
   enableforms;
end;


procedure independentRun;
// this function is called if no file Project1.ini exists
// runs parameter scans by reading in ini files from 1.ini to 100.ini
var i : integer;

begin
for i := 1 to 100 do
begin
  randomize;
  inifilename := ExtractFilePath(ParamStr(0))+ inttostr(i)+'.ini';
  if FileExists(inifilename)then
  begin
    read_ini;
    // security settings
    showvalues := false;
    showwarnings := false;
    calculations := false;
    // initialization
    randomize_patches(0);
    randomize_costs(0);
    statisticsfilename := statisticsfilename + '_Lauf_'+ inttostr(i) ;
    case variationchoice of
      0:
      begin
          if optimizationon then  einsdscanoptimization(0,0) else einsdscan(0,0);
      end;
      1:
      begin
          zweidscan(StartDistribution,costvariation);
      end;
      2:
      begin
          optincentivescan(StartDistribution,costvariation);
      end;
      3:
      begin
          einsdscaneco(StartDistribution,costvariation);
      end;
    end;
  end;
end;

// terminate
application.terminate;
end;


// main procedure
begin
randomize;
forcedirectories(ExtractFilePath(ParamStr(0))+'\statistics');

inifilename := ExtractFilePath(ParamStr(0))+ 'Project1.ini';
if (not FileExists(inifilename) ) then independentrun;

read_ini;
randomize_patches(0);
randomize_costs(0);



end.






