program Project1;

uses
  Forms,
  Unit1mitpreis in 'Unit1mitpreis.pas' {Form1};

{$E exe}

{$R *.RES}

begin
  Application.Initialize;
  Application.Title := 'LatticeTrade';
  Application.CreateForm(TForm1, Form1);
  Form1.Button1Click(form1);
  Application.Run;

end.
