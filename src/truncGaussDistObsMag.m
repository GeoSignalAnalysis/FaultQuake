function [dist_obs_mag]=truncGaussDistObsMag(pdf_magnitudes, x_range_of_mag, mag,sdmag)

dist_mag=pdf_magnitudes(end,:);
dist_obs_mag(size(dist_mag))=zeros;
lower_threshold=(mag-sdmag);
upper_threshold=(mag+sdmag);
lowest_value=find(x_range_of_mag <= lower_threshold,1,'last');
highest_value=find(x_range_of_mag <= upper_threshold,1,'last');

dist_obs_mag(1,lowest_value:highest_value)=dist_mag(1,lowest_value:highest_value);