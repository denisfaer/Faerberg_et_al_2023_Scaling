program wt_permutate;


const
  k = 10000;
  n = 125;
  troubleshoot = FALSE;
  fact = 0.5;
  dau_min = 700;
  dau_max = 1150;
  adult_min = 1900;
  adult_max = 2550;

var
  fin, fout1, fout2, fout3, fout4, fout5: text;
  j, i, y, i1, i2, i3, i4, track, tau: integer;
  s, tmp: string;
  count1, count2, maxT, minT, avg, stdev, cov, temp: real;
  val: array[1..n, 1..5]of real;
  back, rand: array[1..n, 1..5]of real;
  traj: array[1..k, 1..5]of real;


begin
  Randomize;
  assign(fout1, 'dauer_surv.csv');
  rewrite(fout1);
  assign(fout2, 'dauer_curve.csv');
  rewrite(fout2);
  assign(fout3, 'adult_surv.csv');
  rewrite(fout3);
  assign(fout4, 'adult_curve.csv');
  rewrite(fout4);
  assign(fout5, 'normality.csv');
  rewrite(fout5);
  
  //import empirical stages
  assign(fin, 'stages.csv');
  reset(fin);
  readln(fin, s);
  for i := 1 to n do
  begin
    readln(fin, s);
    track := 1;
    tmp := '';
    for y := 1 to length(s) do
      if(s[y] <> ',') then
        tmp := tmp + s[y]
      else
      begin
        back[i, track] := StrToReal(tmp);
        tmp := '';
        inc(track);
      end;
    back[i, track] := StrToReal(tmp);
    back[i, 5] := back[i, 1] + back[i, 2] + back[i, 3] + back[i, 4];
  end;
  close(fin);
  
  //generate randomized trajectory totals
  for i := 1 to k do
  begin
    temp := back[random(n) + 1, 1] + back[random(n) + 1, 2] + back[random(n) + 1, 3] + back[random(n) + 1, 4];
    writeln(fout5, temp);
  end;
  close(fout5);
  
  // dauer commitment exit section
  
  writeln(fout1, '"MinT","MaxT","Window"');
  for y := 1 to k do
  begin
    // generate a random plate
    for j := 1 to n do
    begin
      for i := 1 to 4 do
        rand[j, i] := back[random(n) + 1, i];
      rand[j, 5] := rand[j, 1] + rand[j, 2] + rand[j, 3] + rand[j, 4];
    end;
    
    // maxT and minT
    maxT := 0;
    for i := 1 to n do
      if((rand[i, 1] + fact * rand[i, 2]) > maxT) then
        maxT := rand[i, 1] + fact * rand[i, 2];
    minT := 99999;
    for i := 1 to n do
      if((rand[i, 1] + fact * rand[i, 2]) < minT) then
        minT := rand[i, 1] + fact * rand[i, 2];
    
    writeln(fout1, minT, ',', maxT, ',', maxT - minT);
    
    // generate a sample pair of curves
    if(y = 1) then
    begin
      writeln(fout2, '"Minute","Empyrical","Randomized"');
      for j := dau_min to dau_max do
      begin
        count1 := 0;
        for i := 1 to n do
          if((back[i, 1] + back[i, 2] * fact) > j) then
            count1 := count1 + 1;
        count1 := count1 * 100 / n;
        count2 := 0;
        for i := 1 to n do
          if((rand[i, 1] + rand[i, 2] * fact) > j) then
            count2 := count2 + 1;
        count2 := count2 * 100 / n;
        writeln(fout2, j, ',', count1, ',', count2);
      end;
      close(fout2);
    end;
  end;
  close(fout1);
  
  // maxT and minT of empirical set [dauer]
  maxT := 0;
  for i := 1 to n do
    if((back[i, 1] + fact * back[i, 2]) > maxT) then
      maxT := back[i, 1] + fact * back[i, 2];
  writeln('Dauer MaxT = ', maxT);
  minT := 99999;
  for i := 1 to n do
    if((back[i, 1] + fact * back[i, 2]) < minT) then
      minT := back[i, 1] + fact * back[i, 2];
  writeln('Dauer MinT = ', minT);
  writeln('Dauer Window = ', maxT - minT);
  
  
  // L4/Adult section
  
  writeln(fout3, '"MinT","MaxT","Window"');
  for y := 1 to k do
  begin
    // generate a random plate
    for j := 1 to n do
    begin
      for i := 1 to 4 do
        rand[j, i] := back[random(n) + 1, i];
      rand[j, 5] := rand[j, 1] + rand[j, 2] + rand[j, 3] + rand[j, 4];
    end;
    
    // maxT and minT
    maxT := 0;
    for i := 1 to n do
      if(rand[i, 5] > maxT) then
        maxT := rand[i, 5];
    minT := 99999;
    for i := 1 to n do
      if(rand[i, 5] < minT) then
        minT := rand[i, 5];
    
    writeln(fout3, minT, ',', maxT, ',', maxT - minT);
    
    // generate a sample pair of curves
    if(y = 1) then
    begin
      writeln(fout4, '"Minute","Empyrical","Randomized"');
      for j := adult_min to adult_max do
      begin
        count1 := 0;
        for i := 1 to n do
          if(back[i, 5] < j) then
            count1 := count1 + 1;
        count1 := count1 * 100 / n;
        count2 := 0;
        for i := 1 to n do
          if(rand[i, 5] < j) then
            count2 := count2 + 1;
        count2 := count2 * 100 / n;
        writeln(fout4, j, ',', count1, ',', count2);
      end;
      close(fout4);
    end;
  end;
  close(fout3);
  
  // maxT and minT of empirical set [adult]
  maxT := 0;
  for i := 1 to n do
    if(back[i, 5] > maxT) then
      maxT := back[i, 5];
  writeln('Adult MaxT = ', maxT);
  minT := 99999;
  for i := 1 to n do
    if(back[i, 5] < minT) then
      minT := back[i, 5];
  writeln('Adult MinT = ', minT);
  writeln('Adult Window = ', maxT - minT);
end.