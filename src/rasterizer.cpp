#include "rasterizer.h"

using namespace std;

namespace CGL {

RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
                             size_t width, size_t height,
                             unsigned int sample_rate) {
  this->psm = psm;
  this->lsm = lsm;
  this->width = width;
  this->height = height;
  this->sample_rate = sample_rate;

  sample_buffer.resize(width * height * sample_rate, Color::White);
}

// Used by rasterize_point and rasterize_line
void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
  // TODO: Task 2: You might need to this function to fix points and lines (such
  // as the black rectangle border in test4.svg) NOTE: You are not required to
  // implement proper supersampling for points and lines It is sufficient to use
  // the same color for all supersamples of a pixel for points and lines (not
  // triangles)

  int n = (int)sqrt(sample_rate);
  for (int a = 0; a < n; a++) {
    for (int b = 0; b < n; b++) {
      int index = (y * width + x) * sample_rate + (a * n + b);
      sample_buffer[index] = c;
    }
  }
  // sample_buffer[y * width + x] = c;
}

// Rasterize a point: simple example to help you start familiarizing
// yourself with the starter code.
//
void RasterizerImp::rasterize_point(float x, float y, Color color) {
  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width)
    return;
  if (sy < 0 || sy >= height)
    return;

  fill_pixel(sx, sy, color);
  return;
}

// Rasterize a line.
void RasterizerImp::rasterize_line(float x0, float y0, float x1, float y1,
                                   Color color) {
  if (x0 > x1) {
    swap(x0, x1);
    swap(y0, y1);
  }

  float pt[] = {x0, y0};
  float m = (y1 - y0) / (x1 - x0);
  float dpt[] = {1, m};
  int steep = abs(m) > 1;
  if (steep) {
    dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
    dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
  }

  while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
    rasterize_point(pt[0], pt[1], color);
    pt[0] += dpt[0];
    pt[1] += dpt[1];
  }
}

// Rasterize a triangle.
void RasterizerImp::rasterize_triangle(float x0, float y0, float x1, float y1,
                                       float x2, float y2, Color color) {
  // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
  // valid(0, 0);
  // cout << "sample rate: " << sample_rate << "\n";
  // cout << "height: " << height << " width: " << width << "\n";
  int min_x = (int)floor(min({x0, x1, x2}));
  int max_x = (int)ceil(max({x0, x1, x2}));
  int min_y = (int)floor(min({y0, y1, y2}));
  int max_y = (int)ceil(max({y0, y1, y2}));
  // cout << " minx: " << min_x << " max_x: " << max_x << "  min_y: " << min_y<<
  // " max_y: " << max_y << "\n";
  for (int i = min_x; i <= max_x; i++) {
    for (int j = min_y; j <= max_y; j++) {
      int n = sqrt(sample_rate);
      // Offset is 1 / (2n)
      // Add half to the points you care about
      for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
          float offset = 1.0f / (2 * n);
          float start_x = (a * (1.0f / n));
          float start_y = (b * (1.0f / n));
          float x = start_x + offset + i;
          float y = start_y + offset + j;
          int index = (j * width + i) * sample_rate + (a * n + b);
          auto valid = [&](float x0, float y0, float x1, float y1) {
            float res = (x - x0) * -(y1 - y0);
            res += (y - y0) * (x1 - x0);
            return res >= 0;
          };
          bool clockwise = 1, counter_clockwise = 1;
          clockwise &= valid(x0, y0, x1, y1);
          clockwise &= valid(x1, y1, x2, y2);
          clockwise &= valid(x2, y2, x0, y0);

          counter_clockwise &= valid(x0, y0, x2, y2);
          counter_clockwise &= valid(x2, y2, x1, y1);
          counter_clockwise &= valid(x1, y1, x0, y0);
          if (clockwise || counter_clockwise) {
            sample_buffer[index] = color;
          }
        }
      }
      Color sum_color(0, 0, 0);
      for (int s = 0; s < sample_rate; s++) {
        sum_color += sample_buffer[(j * width + i) * sample_rate + s];
      }
      fill_pixel(i, j, sum_color * (1.0f / sample_rate));
    }
  }
  // TODO: Task 2: Update to implement super-sampled rasterization
}

void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0,
                                                          Color c0, float x1,
                                                          float y1, Color c1,
                                                          float x2, float y2,
                                                          Color c2) {
  // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates
  // and using them to interpolate vertex colors across the triangle Hint: You
  // can reuse code from rasterize_triangle
  int min_x = (int)floor(min({x0, x1, x2}));
  int max_x = (int)ceil(max({x0, x1, x2}));
  int min_y = (int)floor(min({y0, y1, y2}));
  int max_y = (int)ceil(max({y0, y1, y2}));
  // cout << " minx: " << min_x << " max_x: " << max_x << "  min_y: " << min_y<<
  // " max_y: " << max_y << "\n";
  for (int x = min_x; x <= max_x; x++) {
    for (int y = min_y; y <= max_y; y++) {
      float d = -(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1);
      float alpha = -(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1);
      alpha /= d;
      float beta = -(x - x2) * (y0 - y2) + (y - y2) * (x0 - x2);
      beta /= d;
      float gamma = 1.0f - alpha - beta;

      if (alpha >= 0 && beta >= 0 && gamma >= 0) {
        fill_pixel(x, y, alpha * c0 + beta * c1 + gamma * c2);
      }
    }
  }
}

void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0,
                                                float v0, float x1, float y1,
                                                float u1, float v1, float x2,
                                                float y2, float u2, float v2,
                                                Texture &tex) {
  // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample
  // function.
  float min_x = floor(std::min({x0, x1, x2}));
  float max_x = ceil(std::max({x0, x1, x2}));
  float min_y = floor(std::min({y0, y1, y2}));
  float max_y = ceil(std::max({y0, y1, y2}));

  float denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
  for (int x = (int)min_x; x < (int)max_x; x++) {
    for (int y = (int)min_y; y < (int)max_y; y++) {
      float d = -(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1);
      float alpha = -(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1);
      alpha /= d;
      float beta = -(x - x2) * (y0 - y2) + (y - y2) * (x0 - x2);
      beta /= d;
      float gamma = 1.0f - alpha - beta;

      if (alpha >= 0 && beta >= 0 && gamma >= 0) {
        float u = alpha * u0 + beta * u1 + gamma * u2;
        float v = alpha * v0 + beta * v1 + gamma * v2;
        SampleParams sp;
        sp.p_uv = Vector2D(u, v);
        Color color;
        if (psm == P_NEAREST) {
          color = tex.sample_nearest(sp.p_uv, 0);
        } else {
          color = tex.sample_bilinear(sp.p_uv, 0);
        }
        fill_pixel(x, y, color);
      }
    }
  }
  // TODO: Task 6: Set the correct barycentric differentials in the SampleParams
  // struct. Hint: You can reuse code from
  // rasterize_triangle/rasterize_interpolated_color_triangle
}

void RasterizerImp::set_sample_rate(unsigned int rate) {
  // TODO: Task 2: You may want to update this function for supersampling
  // support

  this->sample_rate = rate;
  clear_buffers();
  this->sample_buffer.resize(width * height * sample_rate, Color::White);
}

void RasterizerImp::set_framebuffer_target(unsigned char *rgb_framebuffer,
                                           size_t width, size_t height) {
  // TODO: Task 2: You may want to update this function for supersampling
  // support

  this->width = width;
  this->height = height;
  this->rgb_framebuffer_target = rgb_framebuffer;

  this->sample_buffer.resize(width * height * sample_rate, Color::White);
}

void RasterizerImp::clear_buffers() {
  std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height,
            255);
  std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
}

// This function is called at the end of rasterizing all elements of the
// SVG file.  If you use a supersample buffer to rasterize SVG elements
// for antialising, you could use this call to fill the target framebuffer
// pixels from the supersample buffer data.
//
void RasterizerImp::resolve_to_framebuffer() {
  // TODO: Task 2: You will likely want to update this function for
  // supersampling support
  for (int x = 0; x < width; x++) {
    for (int y = 0; y < height; y++) {
      Color sum_color(0, 0, 0);
      for (int i = 0; i < sample_rate; i++) {
        sum_color += sample_buffer[(y * width + x) * sample_rate + i];
      }
      Color c = sum_color * (1.0f / sample_rate);

      int idx = 3 * (y * width + x);
      rgb_framebuffer_target[idx] = c.r * 255;
      rgb_framebuffer_target[idx + 1] = c.g * 255;
      rgb_framebuffer_target[idx + 2] = c.b * 255;
    }
  }
}

Rasterizer::~Rasterizer() {}

} // namespace CGL
