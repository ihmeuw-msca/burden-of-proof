library(pdftools)
library(magick)
library(gridExtra)
library(grid)
library(png)

list_of_risk_acause_pairs<- c("low_hgb-_all", 
                              "low_hgb-maternal_hem")
list_of_title_pairs<- c("Low Hemoglobin - Maternal Mortality",
                        "Low Hemoglobin - Maternal Hemorrhage")

for (i in list_of_risk_acause_pairs){
  filepath<- paste0("PATH TO PARENT FOLDER",
                    list_of_risk_acause_pairs[i],
                    "/",
                    list_of_risk_acause_pairs[i],
                    ".pdf")
  plot<- magick::image_read_pdf(filepath, density = 300)
  plot<- magick::image_border(plot, "white", "50x100")
  plot<- magick::image_annotate(
    plot,
    text = paste0("Relative Risk Curve for the",  list_of_title_pairs[i],
    size = 8,
    gravity = "north",
    location = "+50+50",
    color = "black"
  )
}

felepath<- paste0("PATH TO PARENT FOLDER",
                  "low_hgb-_all",
                  "/",
                  "low_hgb-_all",
                  ".pdf")

p0 <- magick::image_read_pdf('PATH TO PDF', density = 300)
p0
img_with_margins <- image_border(p0, "white", "50x100")
img_with_margins
image_with_title <- image_annotate(
  img_with_margins,
  text = "Relative Risk Curve for the Low Hemoglobin - Peripartum Depression",
  size = 8,
  gravity = "north",
  location = "+50+50",
  color = "black"
)
image_with_title

image_with_title2 <- image_annotate(
  img_with_margins,
  text = "Relative Risk Curve for the Low Hemoglobin - Maternal Hemorrhage",
  size = 8,
  gravity = "north",
  location = "+50+50",
  color = "black"
)

image_with_title3 <- image_annotate(
  img_with_margins,
  text = "Relative Risk Curve for the Low Hemoglobin - SS",
  size = 8,
  gravity = "north",
  location = "+50+50",
  color = "black"
)
image_with_title3
combined_image <- image_append(c(image_with_title, image_with_title2, image_with_title3), stack = TRUE)
combined_image
test_combined<-grid.arrange(
                    rasterGrob(image_with_title),
                   rasterGrob(image_with_title2),
                   rasterGrob(image_with_title3),
                   ncol=2)
test_combined
ggsave(file="PATH TO PDF", test_combined) #saves g
