function ok = fetch_outcome_data(opts)
% FETCH_OUTCOME_DATA  Download (or verify) the four monthly outcome series.
%
% TARGET FILES (all in data/raw/, SDMX-CSV with TIME_PERIOD + OBS_VALUE)
% ----------------------------------------------------------------------
%   ip_ea.csv      euro-area industrial production, volume index,
%                  industry excl. construction (B-D), seasonally and
%                  calendar adjusted.        [Eurostat sts_inpr_m]
%   hicp_ea.csv    euro-area HICP all-items index (2015 = 100).
%                                             [Eurostat prc_hicp_midx]
%   rate1y_ea.csv  1-year nominal rate LEVEL, monthly average.  Baseline:
%                  12-month EURIBOR (available from 1999; caveat: contains
%                  a bank credit premium, material 2008-2012 -- documented
%                  in docs/CH7_DESIGN.md Sec. 2).  Preferred if you can
%                  locate it in the ECB portal: a 1-year OIS/EONIA-swap
%                  level series.               [ECB Data Portal, FM.*]
%   stoxx50_ea.csv EURO STOXX 50 price index, monthly close/average level.
%                                              [ECB Data Portal, FM.*]
%
% BEHAVIOUR
% ---------
% For each series, if the file already exists it is validated and kept.
% Otherwise the candidate URLs below are tried in order (websave, falling
% back to urlwrite for older Octave).  If all downloads fail -- offline
% machine, changed series keys -- MANUAL INSTRUCTIONS are printed; the rest
% of the pipeline only requires that the four files exist and parse.
%
% IMPORTANT: the series keys below were written from documentation and are
% flagged "verify": portals occasionally re-key series.  If a URL 404s,
% follow the printed manual route (portal search terms given), export as
% SDMX-CSV / csvdata, and save under the target filename.  read_sdmx_csv
% only needs TIME_PERIOD and OBS_VALUE columns, whatever else is present.

if nargin < 1, opts = struct(); end
if ~isfield(opts, 'raw_dir'), opts.raw_dir = fullfile('data', 'raw'); end
if ~exist(opts.raw_dir, 'dir'), mkdir(opts.raw_dir); end

eust = 'https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/';
ecb  = 'https://data-api.ecb.europa.eu/service/data/';

spec = struct('file', {}, 'urls', {}, 'manual', {});

spec(1).file = 'ip_ea.csv';
spec(1).urls = { ...
  [eust 'sts_inpr_m/M.PROD.B-D.SCA.I15.EA20/?format=SDMX-CSV&startPeriod=1999-01'], ...
  [eust 'sts_inpr_m/M.PROD.B-D.SCA.I15.EA19/?format=SDMX-CSV&startPeriod=1999-01']};
spec(1).manual = ['Eurostat Data Browser, dataset sts_inpr_m: Production in ' ...
  'industry, monthly.  Filter: NACE B-D (industry excl construction), ' ...
  'seasonally and calendar adjusted (SCA), index 2015=100 (I15), geo EA20 ' ...
  '(or EA19), from 1999-01.  Download as SDMX-CSV.'];

spec(2).file = 'hicp_ea.csv';
spec(2).urls = { ...
  [eust 'prc_hicp_midx/M.I15.CP00.EA/?format=SDMX-CSV&startPeriod=1999-01'], ...
  [eust 'prc_hicp_midx/M.I15.CP00.EA19/?format=SDMX-CSV&startPeriod=1999-01']};
spec(2).manual = ['Eurostat Data Browser, dataset prc_hicp_midx: HICP monthly ' ...
  'index.  Filter: all-items CP00, index 2015=100 (I15), geo EA (changing ' ...
  'composition), from 1999-01.  Download as SDMX-CSV.  (NSA index; see ' ...
  'CH7_DESIGN Sec. 2 for the seasonality discussion.)'];

spec(3).file = 'rate1y_ea.csv';
spec(3).urls = { ...
  [ecb 'FM/M.U2.EUR.RT.MM.EURIBOR1YD_.HSTA?format=csvdata&startPeriod=1999-01'], ...
  [ecb 'FM/M.U2.EUR.RT.MM.EURIBOR1YD.HSTA?format=csvdata&startPeriod=1999-01']};
spec(3).manual = ['ECB Data Portal (data.ecb.europa.eu), search "EURIBOR ' ...
  '12-month  monthly": Euribor 1-year, historical close, average of ' ...
  'observations through period, euro area.  Export CSV (csvdata).  If you ' ...
  'find a 1-year OIS / EONIA swap LEVEL series with 1999+ coverage, use it ' ...
  'instead under the same filename and note the swap in the thesis.'];

spec(4).file = 'stoxx50_ea.csv';
spec(4).urls = { ...
  [ecb 'FM/M.U2.EUR.DS.EI.DJES50I.HSTA?format=csvdata&startPeriod=1999-01']};
spec(4).manual = ['ECB Data Portal, search "Dow Jones EURO STOXX 50  monthly": ' ...
  'price index, historical close, average of observations.  Export CSV.'];

ok = true;
for i = 1:numel(spec)
    path = fullfile(opts.raw_dir, spec(i).file);
    if exist(path, 'file')
        try
            [ym, ~] = read_sdmx_csv(path);
            fprintf('fetch_outcome_data: %-14s present, %d monthly obs (%d-%02d to %d-%02d)\n', ...
                    spec(i).file, numel(ym), floor((ym(1)-1)/12), ym(1)-12*floor((ym(1)-1)/12), ...
                    floor((ym(end)-1)/12), ym(end)-12*floor((ym(end)-1)/12));
            continue
        catch
            fprintf('fetch_outcome_data: %s exists but does not parse; re-downloading.\n', spec(i).file);
        end
    end
    got = false;
    for u = 1:numel(spec(i).urls)
        try
            if exist('websave', 'file') || exist('websave', 'builtin')
                websave(path, spec(i).urls{u});
            else
                urlwrite(spec(i).urls{u}, path);  %#ok<URLWR>
            end
            [~, ~] = read_sdmx_csv(path);           % validate
            fprintf('fetch_outcome_data: downloaded %s\n', spec(i).file);
            got = true;
            break
        catch
            if exist(path, 'file'), delete(path); end
        end
    end
    if ~got
        ok = false;
        fprintf(2, ['\nfetch_outcome_data: could NOT obtain %s automatically.\n' ...
                    'MANUAL ROUTE: %s\nSave the file as %s\n\n'], ...
                spec(i).file, spec(i).manual, path);
    end
end
if ok
    fprintf('fetch_outcome_data: all four outcome files ready.\n');
else
    fprintf(2, 'fetch_outcome_data: some files missing -- follow the manual instructions above, then rerun.\n');
end
end
