unit SAOptimization;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  TeEngine, Series, ExtCtrls, TeeProcs, Chart;

type
  TSimmulatedAnnealing = class(TForm)
    Chart1: TChart;
    Series1: TFastLineSeries;
    Chart2: TChart;
    Series2: TFastLineSeries;
    Chart3: TChart;
    Series3: TFastLineSeries;
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  SimmulatedAnnealing: TSimmulatedAnnealing;

implementation

uses Unit1mitpreis;

{$R *.DFM}




end.
