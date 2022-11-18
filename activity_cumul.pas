program wt_l3_activity;


const
  k = 125;
  red = 30;
  fps = 3;
  name = 'disp.csv';
  maxval = 30000;
  summary = TRUE;
  bins = 10;    //number of large bins
  border = 1;   //hours to ignore on L3 sides
  binrep = 10000;   //number of small bins runs repetitions
  bincons = 10; //number of consucutive small bins
  smallbin = 5;//minutes in a small bin

var
  fdisp, fout, fin: text;
  j, i, t, track, bin, y: integer;
  ts, move, disp, avg, l1, l2, l3: real;
  val, ds: array [1..k, 1..maxval]of real;
  max: array[1..k]of integer;
  range: array[1..k, 1..2]of integer;
  onerun: array[1..k, 1..bincons]of real;
  allruns: array[1..binrep, 2..bincons]of real;
  aver, stdev: array [1..bincons]of real;
  s: string;

function correl(x, y: integer): real;
var
  cnt: real;
  u: integer;
begin
  cnt := 0;
  for u := 1 to k do
    cnt := cnt + (onerun[u, x] - aver[x]) * (onerun[u, y] - aver[y]) / (stdev[x] * stdev[y]);
  correl := cnt / (k - 1);
end;

begin
  //reset the arrays
  for j := 1 to k do
    for i := 1 to maxval do
      val[j, i] := -1;
  for j := 1 to k do
    for i := 1 to maxval do
      ds[j, i] := -1;
  for j := 1 to k do
    max[j] := 0;
  //generate cumulative curves
  for j := 1 to k do
  begin
    move := 0;
    track := 0;
    assign(fdisp, IntToStr(j) + 're' + IntToStr(red) + name);
    reset(fdisp);
    while(not eof(fdisp)) do
    begin
      readln(fdisp, disp);
      move := move + disp;
      inc(track);
      val[j, track] := move;
      ds[j, track] := disp;
    end;
    max[j] := track;
    close(fdisp);
  end;
  //read L3 boundaries
  assign(fin, 'stages.csv');
  reset(fin);
  readln(fin, s);
  for i := 1 to k do
  begin
    readln(fin, s); 
    l1 := StrToFloat(copy(s, 1, pos(',', s) - 1));
    delete(s, 1, pos(',', s));
    l2 := StrToFloat(copy(s, 1, pos(',', s) - 1));
    delete(s, 1, pos(',', s));
    l3 := StrToFloat(copy(s, 1, pos(',', s) - 1));
    delete(s, 1, pos(',', s));
    range[i, 1] := round((l1 + l2 + border * 60) * fps * 60 / red);
    range[i, 2] := range[i, 1] + round((l3 - border * 120) * fps * 60 / red);
  end;
  close(fin);
  //general summary
  if(summary) then
  begin
    avg := 0;
    for j := 1 to k do
      avg := avg + max[j];
    avg := avg * red / (k * fps * 3600);
    writeln('Average recording: ', avg:3, ' h');
    avg := 999999;
    for j := 1 to k do
      if((max[j] * red / (fps * 3600)) < avg) then
        avg := max[j] * red / (fps * 3600);
    writeln('Shortest recording: ', avg:3, ' h');
    avg := 0;
    for j := 1 to k do
      if((max[j] * red / (fps * 3600)) < 50) then
        avg := avg + 1;
    writeln('Recordings under 50h: ', avg);
    writeln();
    avg := 0;
    for j := 1 to k do
      avg := avg + range[j, 2] - range[j, 1];
    avg := avg * red / (k * fps * 3600);
    writeln('Average L3: ', avg:3, ' h');
    avg := 99999;
    for j := 1 to k do
      if((range[j, 2] - range[j, 1]) < avg) then
        avg := range[j, 2] - range[j, 1];
    avg := avg * red / (fps * 3600); 
    writeln('Shortest L3: ', avg:3, ' h');
    avg := 0;
    for j := 1 to k do
      if((range[j, 2] - range[j, 1]) * red / (fps * 3600) < 5) then
        avg := avg + 1;
    writeln('L3 under 5h: ', avg);
  end;
  //bin L3 into bins
  assign(fout, 'activity_cumul.csv');
  rewrite(fout);
  for j := 1 to k do
    if((range[j, 2] - range[j, 1]) * red / (fps * 3600) >= 5) then
    begin
      bin := round(30 * 60 * fps / red);
      for i := 1 to bins - 1 do
        write(fout, val[j, range[j, 1] + i * bin - 1] - val[j, range[j, 1] + (i - 1) * bin], ',');
      writeln(fout, val[j, range[j, 1] + bins * bin - 1] - val[j, range[j, 1] + (bins - 1) * bin]);
    end;
  close(fout);
  //small bins
  assign(fout, 'smallbins.csv');
  rewrite(fout);
  writeln(fout,'"Correlation","Distance"');
  for y := 1 to binrep do
  begin
    //generate random bin runs
    for j := 1 to k do
    begin
      bin := random(trunc((range[j, 2] - range[j, 1]))) + range[j, 1];
      for i := 1 to bincons do
        onerun[j, i] := val[j, bin + i * round(smallbin * 60 * fps / red) - 1] - val[j, bin + (i - 1) * round(smallbin * 60 * fps / red)];
    end;
    //generate averages and stdevs
    for i := 1 to bincons do
    begin
      aver[i] := 0;
      stdev[i] := 0;
      for j := 1 to k do
        aver[i] := aver[i] + onerun[j, i];
      aver[i] := aver[i] / k;
      for j := 1 to k do
        stdev[i] := stdev[i] + (aver[i] - onerun[j, i]) * (aver[i] - onerun[j, i]);
      stdev[i] := stdev[i] / (k - 1);
      stdev[i] := sqrt(stdev[i]);
    end;
    //save correl values
    for i := 2 to bincons do
      allruns[y, i] := correl(1, i);
    for i:= 2 to bincons do
    writeln(fout,allruns[y,i],',"D',i-1,'"');
   // writeln(fout,allruns[y,bincons]);
  end;
  close(fout);
end.