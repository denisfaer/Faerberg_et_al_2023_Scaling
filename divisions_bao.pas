program embryo_divisions;

const
  n = 50; //number of embryos
  cens = 44;//censor cells where the amount of embryos with observed division is less or equal to this threshold
  cellmax = 1000;
  branchmax = 2000;

var
  names: array [1..cellmax] of string;
  division: array [1..n, 1..cellmax + 1] of real;
  labels: array[1..branchmax] of string;
  branch: array [1..n, 1..branchmax] of real;
  recording: array [1..n] of real;
  gen: array [1..branchmax] of integer;
  perm: array[1..cellmax] of boolean;
  counter, embryo, loc, i, j, confirm, cnt, loc2, bcounter, gen_c, y: integer;
  fin, fout, fout2, fout3, fout4, fout5: text;
  s, temp, cell, lab, nm, trajectory: string;
  divis, cc_temp: real;
  avg: array[1..6]of real;

function find1(x: string): integer;
var
  i_tmp: integer;
begin
  i_tmp := 1;
  while((x <> names[i_tmp]) and (i_tmp < cellmax)) do
    inc(i_tmp);
  result := i_tmp;
end;

function find2(x: string): integer;
var
  i_tmp: integer;
begin
  i_tmp := 1;
  while((x <> labels[i_tmp]) and (i_tmp < branchmax)) do
    inc(i_tmp);
  result := i_tmp;
end;


procedure traj(embr: integer; var x: string);
var
  i_tmp: integer;
  x_tmp: string;
begin
  if(length(x) > 3) then
  begin
    i_tmp := 1;
    while((x <> copy(labels[i_tmp], pos('-', labels[i_tmp]) + 1, length(labels[i_tmp]) - pos('-', labels[i_tmp]))) and (i_tmp < branchmax)) do
      inc(i_tmp);
    x_tmp := copy(x, 1, length(x) - 1);
    trajectory := FloatToStr(branch[embr, i_tmp]) + ',' + trajectory;
    traj(embr, x_tmp);
  end;
end;


// grab cell name from row: copy(s,pos(',',s)+1,pos(',',s,pos(',',s)+1)-pos(',',s)-1)

begin
  assign(fout, 'division_times.csv');
  rewrite(fout);
  assign(fout2, 'branches.csv');
  rewrite(fout2);
  assign(fout3, 'cellcycle_series.csv');
  rewrite(fout3);
  assign(fout4, 'cellcycle_batching.csv');
  rewrite(fout4);
  assign(fout5, 'bao_norm.csv');
  rewrite(fout5);
  for i := 1 to n do
    for j := 1 to cellmax do
      division[i, j] := -1;
  for i := 1 to cellmax do
    names[i] := '';
  for i := 1 to n do
    for j := 1 to branchmax do
      branch[i, j] := -1;
  for i := 1 to branchmax do
    labels[i] := '';
  for i := 1 to branchmax do
    gen[i] := 0;
  for i := 1 to n do
    recording[i] := 0;
  counter := 0;
  bcounter := 0;
  
  for y := 1 to n do
  begin
    assign(fin, 'embryo (' + y + ').csv');
    reset(fin);
    while(not eof(fin)) do
    begin
      readln(fin, s);
      cell := copy(s, 1, pos(',', s) - 1);                                    //grab cell name
      temp := copy(s, pos(',', s) + 1, length(s) - pos(',', s));
      divis := StrToInt(copy(temp, pos(',', temp) + 1, length(temp) - pos(',', temp)));
      loc := find1(cell);                                                     //find index of the cell name in NAMES array
      if(loc < cellmax) then
        division[y, loc] := divis
        else
      begin
        inc(counter);                                                         //create new entry in NAMES array on first cell name encounter
        names[counter] := cell;
        division[y, counter] := divis;
      end;
    end;
    close(fin);
  end;
  
  // censor cells below threshold
  for i := 1 to counter do
  begin
    cnt := 0;
    for j := 1 to n do
      if(division[j, i] <> -1) then
        inc(cnt);
    perm[i] := (cnt > cens);
  end;
  
  // calculate branch lengths
  for i := 1 to n do
    for j := 1 to counter do
      if(division[i, j] > 0) then
        if((perm[j]) and (names[j] <> 'ABa') and (names[j] <> 'ABp') and (names[j] <> 'EMS') and (names[j] <> 'P2')) then
        begin
          cell := names[j];
          cell := copy(cell, 1, length(cell) - 1);
          if(names[j] = 'MS') then cell := 'EMS';
          if(names[j] = 'E') then cell := 'EMS';
          if(names[j] = 'C') then cell := 'P2';
          if(names[j] = 'D') then cell := 'P3';
          if(names[j] = 'P3') then cell := 'P2';
          if(names[j] = 'P4') then cell := 'P3';
          if(copy(cell, 1, 2) = 'AB') then gen_c := length(cell) - 2
          else
          if((names[j] = 'MS') or (names[j] = 'E') or (names[j] = 'C') or (names[j] = 'P3')) then gen_c := 1
          else 
          if(copy(names[j], 1, 2) = 'MS') then gen_c := length(names[j]) - 1
          else 
          if(copy(names[j], 1, 1) = 'E') then gen_c := length(names[j])
          else 
          if(copy(names[j], 1, 1) = 'C') then gen_c := length(names[j])
          else 
          if((names[j] = 'D') or (names[j] = 'P4')) then gen_c := 2
          else 
          if(copy(names[j], 1, 1) = 'D') then gen_c := length(names[j]) + 1;
          loc := find1(cell);
          lab := cell + '-' + names[j];
          loc2 := find2(lab);
          if(loc2 < branchmax) then
            branch[i, loc2] :=  division[i, j] - division[i, loc]
          else
          begin
            inc(bcounter);
            labels[bcounter] := lab;
            gen[bcounter] := gen_c;
            branch[i, bcounter] := division[i, j] - division[i, loc];
          end;
        end;
  
  // DIVISION
  // write column labels
  write(fout, '"Embryo",');
  for i := 1 to (counter - 1) do
    if(perm[i]) then
      write(fout, '"', names[i], '",');
  if(perm[counter]) then
    writeln(fout, '"', names[counter], '"')
  else writeln(fout);
  // write division times
  for j := 1 to n do
  begin
    write(fout, j, ',');
    for i := 1 to (counter - 1) do
      if(perm[i]) then
        if(division[j, i] <> -1) then
          write(fout, division[j, i], ',')
        else write(fout, ',');
    if((division[j, counter] <> -1) and (perm[counter])) then
      writeln(fout, division[j, counter])
    else writeln(fout);
  end;
  
  // BRANCHES
  // write column labels
  write(fout2, '"Embryo",');
  for i := 1 to (bcounter - 1) do
    write(fout2, '"', labels[i], '",');
  writeln(fout2, '"', labels[bcounter], '"');
  write(fout2, '"Gen",');
  for i := 1 to (bcounter - 1) do
    write(fout2, '"', gen[i], '",');
  writeln(fout2, '"', gen[bcounter], '"');
  // write branch lengths
  for j := 1 to n do
  begin
    write(fout2, j, ',');
    for i := 1 to (bcounter - 1) do
      if(branch[j, i] > 0) then
        write(fout2, branch[j, i], ',')
      else write(fout2, ',');
    if(branch[j, bcounter] > 0) then
      writeln(fout2, branch[j, bcounter])
    else writeln(fout2);
  end;
  
  // CELL CYCLE SERIES
  // extract branches as genc-part series
  writeln(fout3, '"Embryo","Terminus","CC1","CC2","CC3","CC4"');
  for i := 1 to n do
    for j := 1 to bcounter do
      if((gen[j] = 6) and (copy(labels[j], 1, 2) = 'AB') and (branch[i, j] > 0)) then
      begin
        write(fout3, i, ',', copy(labels[j], pos('-', labels[j]) + 1, length(labels[j]) - pos('-', labels[j])), ',');
        nm := copy(labels[j], 1, pos('-', labels[j]) - 1);
        trajectory := FloatToStr(branch[i, j]);
        traj(i, nm);
        writeln(fout3, trajectory);
      end;
  
  // BATCHING
  // save batching averages
  writeln(fout4, '"Batch","Terminus","Total Average","Fast"');
  for j := 1 to bcounter do
    if((gen[j] = 4) and (copy(labels[j], 1, 2) = 'AB') and (branch[i, j] > 0)) then
    begin
      for y := 1 to 5 do
      begin
        avg[y] := 0;
        //writeln(y, ' - ', copy(labels[j], pos('-', labels[j]) + 1, length(labels[j]) - pos('-', labels[j])));
        for i := (10 * (y - 1) + 1) to 10 * y do
        begin
          nm := copy(labels[j], 1, pos('-', labels[j]) - 1);
          trajectory := FloatToStr(branch[i, j]);
          traj(i, nm);
          //write(trajectory, ' ');
          avg[y] := avg[y] + StrToFloat(copy(trajectory, 1, pos(',', trajectory) - 1));
          delete(trajectory, 1, pos(',', trajectory));
          avg[y] := avg[y] + StrToFloat(copy(trajectory, 1, pos(',', trajectory) - 1));
          delete(trajectory, 1, pos(',', trajectory));
          avg[y] := avg[y] + StrToFloat(copy(trajectory, 1, pos(',', trajectory) - 1));
          delete(trajectory, 1, pos(',', trajectory));
          avg[y] := avg[y] + StrToFloat(trajectory);
          //writeln(avg);
        end;
        avg[y] := avg[y] / 10;
      end;
      avg[6] := (avg[1] + avg[2] + avg[3] + avg[4] + avg[5]) / 5;
      for y := 1 to 5 do
        if(avg[y] < avg[6]) then
          writeln(fout4, y, ',', copy(labels[j], pos('-', labels[j]) + 1, length(labels[j]) - pos('-', labels[j])), ',', avg[y], ',1')
        else writeln(fout4, y, ',', copy(labels[j], pos('-', labels[j]) + 1, length(labels[j]) - pos('-', labels[j])), ',', avg[y], ',0');
    end;
  
  // DE-BATCHING
  // save termini normalized by batch average
  writeln(fout5, '"Batch","Embryo","Terminus","CC1","CC2","CC3","CC4"');
  for y := 1 to 5 do
  begin
    avg[y] := 0;
    counter := 0;
    for i := (10 * (y - 1) + 1) to 10 * y do
      for j := 1 to bcounter do
        if((gen[j] = 4) and (copy(labels[j], 1, 2) = 'AB') and (branch[i, j] > 0)) then
        begin
          nm := copy(labels[j], 1, pos('-', labels[j]) - 1);
          trajectory := FloatToStr(branch[i, j]);
          traj(i, nm);
          avg[y] := avg[y] + StrToFloat(copy(trajectory, 1, pos(',', trajectory) - 1));
          delete(trajectory, 1, pos(',', trajectory));
          avg[y] := avg[y] + StrToFloat(copy(trajectory, 1, pos(',', trajectory) - 1));
          delete(trajectory, 1, pos(',', trajectory));
          avg[y] := avg[y] + StrToFloat(copy(trajectory, 1, pos(',', trajectory) - 1));
          delete(trajectory, 1, pos(',', trajectory));
          avg[y] := avg[y] + StrToFloat(trajectory);
          inc(counter);
        end;
    avg[y] := avg[y] / counter;
    for i := (10 * (y - 1) + 1) to 10 * y do
      for j := 1 to bcounter do
        if((gen[j] = 4) and (copy(labels[j], 1, 2) = 'AB') and (branch[i, j] > 0)) then
        begin
          write(fout5, y, ',', i, ',', copy(labels[j], pos('-', labels[j]) + 1, length(labels[j]) - pos('-', labels[j])), ',');
          nm := copy(labels[j], 1, pos('-', labels[j]) - 1);
          trajectory := FloatToStr(branch[i, j]);
          traj(i, nm);
          cc_temp := StrToFloat(copy(trajectory, 1, pos(',', trajectory) - 1)) / avg[y];
          write(fout5, cc_temp, ',');
          delete(trajectory, 1, pos(',', trajectory));
          cc_temp := StrToFloat(copy(trajectory, 1, pos(',', trajectory) - 1)) / avg[y];
          write(fout5, cc_temp, ',');
          delete(trajectory, 1, pos(',', trajectory));
          cc_temp := StrToFloat(copy(trajectory, 1, pos(',', trajectory) - 1)) / avg[y];
          write(fout5, cc_temp, ',');
          delete(trajectory, 1, pos(',', trajectory));
          cc_temp := StrToFloat(trajectory) / avg[y];
          writeln(fout5, cc_temp);
        end;
  end;
  
  
  close(fout);
  close(fout2);
  close(fout3);
  close(fout4);
  close(fout5);
end.