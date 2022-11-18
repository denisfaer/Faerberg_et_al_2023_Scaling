program mammillata;

const
  n = 10;
  max_cell = 10000;
  min_per_frame = 2;
  max_lin = 1000;

var
  i, j, y, cnt, track, k: integer;
  fin, fout: text;
  mother_id, mother_name, child1, child2: array[1..n, 1..max_cell]of string;
  mother_life: array[1..n, 1..max_cell]of integer;
  embryo_cells: array[1..n]of integer;
  s, outp, last: string;
  lineage: array[1..n, 1..max_lin]of string;

//find lineage
function find2(z: integer; name: string): integer;
var
  x, save: integer;
  found: boolean;
begin
  found := FALSE;
  for x := 1 to max_lin do
    if(lineage[z, x] = name) then
    begin
      found := TRUE;
      save := x;
    end;
  if(found) then
    result := save
  else result := -1;
end;

//find mother cell identitiy
function find(x1: integer; name: string): integer;
var
  y1, save: integer;
  found: boolean;
begin
  found := FALSE;
  for y1 := 1 to embryo_cells[x1] do
    if((child1[x1, y1] = name) or (child2[x1, y1] = name)) then
    begin
      save := y1;
      found := TRUE;
    end;
  if(found) and (mother_name[x1, save] <> '0') and (mother_name[x1, save][2] <> '6') then
    result := save
  else result := -1;
end;

//generate a string of predecessor names
function lineage_names(x2, y2: integer; m_name: string): string;
var
  tmp1: integer;
begin
  tmp1 := find(x2, m_name);
  if(tmp1 > 0) then
    
    result := lineage_names(x2, tmp1, mother_id[x2, tmp1]) + '-' + mother_name[x2, y2]
  else result := mother_name[x2, y2];
end;

//generate a string of predecessor lifetimes
function lineage_durs(x3, y3: integer; mo_name: string): string;
var
  tmp2: integer;
begin
  tmp2 := find(x3, mo_name);
  if(tmp2 > 0) then
    result := lineage_durs(x3, tmp2, mother_id[x3, tmp2]) + ',' + IntToStr(mother_life[x3, y3])
  else result := IntToStr(mother_life[x3, y3]);
end;

begin
  //reset all arrays
  last := 'xxx';
  for j := 1 to n do
    for i := 1 to max_cell do
    begin
      mother_id[j, i] := '0';
      mother_name[j, i] := '0';
      child1[j, i] := '0';
      child2[j, i] := '0';
      mother_life[j, i] := 0;
    end;
  //extract division data from CSV files
  for j := 1 to n do
  begin
    assign(fin, 'embryo (' + IntToStr(j) + ').csv');
    reset(fin);
    readln(fin, s);
    cnt := 0;
    while(not eof(fin)) do
    begin
      readln(fin, s);
      inc(cnt);
      mother_id[j, cnt] := copy(s, 1, pos(',', s) - 1);
      delete(s, 1, pos(',', s));
      delete(s, 1, pos(',', s));
      child1[j, cnt] := copy(s, 1, pos(',', s) - 1);
      delete(s, 1, pos(',', s));
      child2[j, cnt] := copy(s, 1, pos(',', s) - 1);
      delete(s, 1, pos(',', s));
      mother_life[j, cnt] := StrToInt(copy(s, 1, pos(',', s) - 1)) * min_per_frame;
      delete(s, 1, pos(',', s));
      mother_name[j, cnt] := copy(s, 1, pos(',', s) - 1);
    end;
    embryo_cells[j] := cnt;
    close(fin);
  end;
  //generate lineages from terminal cells
  assign(fout, 'mammillata.csv');
  rewrite(fout);
  writeln(fout, '"Embryo","Lineage","CC1","CC2","CC3","CC4"');
  for j := 1 to n do
  begin
    track := 0;
    for i := 1 to embryo_cells[j] do
      if((child1[j, i] = '0') or (child2[j, i] = '0') and (mother_name[j, i] <> '0')) then
      begin
        s := lineage_names(j, i, mother_id[j, i]);
        while((length(s) > 0) and (s[length(s)] <> '-')) do
          delete(s, length(s), 1);
        k := find2(j, s);
        if(k < 0) and (s <> '') then
        begin
          inc(track);
          lineage[j, track] := s;
        end;
        outp := lineage_durs(j, i, mother_id[j, i]);
        while((length(outp) > 0) and (outp[length(outp)] <> ',')) do
          delete(outp, length(outp), 1);
        if(s <> '') and (length(outp) > 4) and (k < 0) then
        begin
          delete(s, length(s), 1);
          writeln(fout, j, ',', s, ',', outp);
        end;
      end;
  end;
  close(fout);
end.