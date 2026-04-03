/// Parse a VoxelDose.csv produced by OutputVoxelDose into a dense 3-D grid.
///
/// CSV format:
///   Row 0: `x_y_z,...` coordinate labels for occupied voxels (f32 µm)
///   Rows 1..N: cumulative dose per image for each occupied voxel
///
/// We take the **last** row as the final dose state.
pub struct Grid {
    pub data: Vec<f32>,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub min: [f32; 3],
    pub spacing: f32,
}

impl Grid {
    pub fn max_dose(&self) -> f32 {
        self.data.iter().cloned().fold(f32::NEG_INFINITY, f32::max)
    }
}

pub fn parse_voxel_dose_csv(content: &str) -> Result<Grid, String> {
    let mut lines = content.lines();

    // --- Header row ---
    let header = lines.next().ok_or("CSV is empty")?;

    let mut coords: Vec<[f32; 3]> = Vec::new();
    for token in header.split(',') {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        let parts: Vec<&str> = token.splitn(3, '_').collect();
        if parts.len() != 3 {
            return Err(format!("unexpected header token: '{token}'"));
        }
        let x: f32 = parts[0]
            .parse()
            .map_err(|_| format!("bad x in '{token}'"))?;
        let y: f32 = parts[1]
            .parse()
            .map_err(|_| format!("bad y in '{token}'"))?;
        let z: f32 = parts[2]
            .parse()
            .map_err(|_| format!("bad z in '{token}'"))?;
        coords.push([x, y, z]);
    }

    if coords.is_empty() {
        return Err("no voxel coordinates in header".into());
    }

    // --- Data rows — keep last ---
    let mut last_row: Vec<f32> = Vec::new();
    for line in lines {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let row: Vec<f32> = line
            .split(',')
            .filter(|s| !s.trim().is_empty())
            .map(|s| s.trim().parse::<f32>().unwrap_or(0.0))
            .collect();
        if !row.is_empty() {
            last_row = row;
        }
    }

    if last_row.is_empty() {
        return Err("no data rows in CSV".into());
    }

    // --- Derive grid dimensions ---
    let xs: Vec<f32> = coords.iter().map(|c| c[0]).collect();
    let ys: Vec<f32> = coords.iter().map(|c| c[1]).collect();
    let zs: Vec<f32> = coords.iter().map(|c| c[2]).collect();

    let min_x = xs.iter().cloned().fold(f32::INFINITY, f32::min);
    let max_x = xs.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let min_y = ys.iter().cloned().fold(f32::INFINITY, f32::min);
    let max_y = ys.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let min_z = zs.iter().cloned().fold(f32::INFINITY, f32::min);
    let max_z = zs.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

    let spacing = min_positive_gap(&xs)
        .or_else(|| min_positive_gap(&ys))
        .or_else(|| min_positive_gap(&zs))
        .unwrap_or(1.0);

    let nx = (((max_x - min_x) / spacing).round() as usize) + 1;
    let ny = (((max_y - min_y) / spacing).round() as usize) + 1;
    let nz = (((max_z - min_z) / spacing).round() as usize) + 1;

    let mut data = vec![0.0f32; nx * ny * nz];

    for (idx, coord) in coords.iter().enumerate() {
        let dose = if idx < last_row.len() {
            last_row[idx]
        } else {
            0.0
        };
        let i = ((coord[0] - min_x) / spacing).round() as usize;
        let j = ((coord[1] - min_y) / spacing).round() as usize;
        let k = ((coord[2] - min_z) / spacing).round() as usize;
        if i < nx && j < ny && k < nz {
            data[i * ny * nz + j * nz + k] = dose;
        }
    }

    Ok(Grid {
        data,
        nx,
        ny,
        nz,
        min: [min_x, min_y, min_z],
        spacing,
    })
}

fn min_positive_gap(vals: &[f32]) -> Option<f32> {
    let mut sorted: Vec<f32> = vals.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    sorted.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
    sorted
        .windows(2)
        .map(|w| w[1] - w[0])
        .filter(|&d| d > 1e-6)
        .fold(None, |acc, d| Some(acc.map_or(d, |prev: f32| prev.min(d))))
}
