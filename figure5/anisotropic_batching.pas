program batching_correl;

const
  steps = 4;
  val_avg = 10;
  val_vari = 1;
  fine = 10;
  scale_max: array[1..fine]of real = (1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5);
  batch_size: array[0..fine]of integer = (1, 2, 5, 10, 20, 35, 50, 75, 100, 200, 500);
  batch_number: array[1..fine - 1]of integer = (2, 3, 4, 5, 6, 7, 8, 9, 10);
  rep = 10000;
  write_sample = FALSE;
  batch = true;

var
  i, j, y, r, p, f1, f2, f3, temp: integer;
  par: real;
  fout, fout1: text;
  scale: array[1..steps] of real;
  val: array[1..10, 1..500, 1..steps]of real;
  avg, stdev: array[1..steps] of real;
  cor: array[1..rep, 1..steps, 1..steps]of real;
  param: array[1..fine, 0..fine, 1..fine - 1, 1..2]of real;


function correl(x, z: integer): real;
var
  cnt: real;
  u1, u2: integer;
begin
  cnt := 0;
  for u1 := 1 to batch_number[f3] do
    for u2 := 1 to batch_size[f2] do
      cnt := cnt + (val[u1, u2, x] - avg[x]) * (val[u1, u2, z] - avg[z]) / (stdev[x] * stdev[z]);
  correl := cnt / (batch_number[f3] * batch_size[f2] - 1);
end;


begin
  randomize();
  //reset output array
  for f1 := 1 to fine do
    for f2 := 0 to fine do
      for f3 := 1 to fine - 1 do
        for y := 1 to 2 do
          param[f1, f2, f3, y] := 0;
  //fill the output array
  for f1 := 1 to fine do
    for f2 := 0 to fine do
      for f3 := 1 to fine - 1 do
      begin
        for r := 1 to rep do
        begin
          //generate randomized sets
          for y := 1 to batch_number[f3] do
          begin
            for i := 1 to steps do
              scale[i] := scale_max[f1] - random() * (scale_max[f1] - 1);
            for j := 1 to batch_size[f2] do
              for i := 1 to steps do
                val[y, j, i] := scale[i] * (val_avg + 2 * (random() - 0.5) * val_vari);
          end;
          //calculate average and stdev for each step
          for i := 1 to steps do
            avg[i] := 0;
          for i := 1 to steps do
            stdev[i] := 0;
          for i := 1 to batch_number[f3] do
            for j := 1 to batch_size[f2] do
              for y := 1 to steps do
                avg[y] := avg[y] + val[i, j, y];
          for i := 1 to steps do
            avg[i] := avg[i] / (batch_number[f3] * batch_size[f2]);
          for i := 1 to batch_number[f3] do
            for j := 1 to batch_size[f2] do
              for y := 1 to steps do
                stdev[y] := stdev[y] + (val[i, j, y] - avg[y]) * (val[i, j, y] - avg[y]);
          for i := 1 to steps do
            stdev[i] := sqrt(stdev[i] / ((batch_number[f3] * batch_size[f2]) - 1));
          //save correl values
          cor[r, 1, 2] := correl(1, 2);
        end;
        //calculate correl average and stdev
        for r := 1 to rep do
          param[f1, f2, f3, 1] := param[f1, f2, f3, 1] + cor[r, 1, 2];
        param[f1, f2, f3, 1] := param[f1, f2, f3, 1] / rep;
        for r := 1 to rep do
          param[f1, f2, f3, 2] := param[f1, f2, f3,  2] + (cor[r, 1, 2] - param[f1, f2, f3,  1]) * (cor[r, 1, 2] - param[f1, f2, f3,  1]);
        param[f1, f2, f3,  2] := sqrt(param[f1, f2, f3,  2] / (rep - 1));     
      end;
  
  //write all repeat correls
  assign(fout1, 'batch_expected.csv');
  rewrite(fout1);
  writeln(fout1, '"Scale Max","Batch Size","Total N","Average","StDev"');
  for f1 := 1 to fine do
    for f2 := 0 to fine do
      for f3 := 1 to fine - 1 do
        writeln(fout1, scale_max[f1], ',', batch_size[f2], ',', batch_size[f2] * batch_number[f3], ',', param[f1, f2, f3,  1], ',', param[f1, f2, f3,  2]);
  close(fout1);
  temp := trunc(Milliseconds / 1000);
  writeln('Run time: ', temp div 60, 'min ', temp mod 60, 'sec');
end.