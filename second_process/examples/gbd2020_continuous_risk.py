from espipeline.main import PostContinuousProcess, ContinuousPipeline
from espipeline.filemanager import FileManager
import warnings
warnings.filterwarnings("ignore")


def main():
    i_folder = "/mnt/team/msca/pub/archive/evidence-score/gbd2020"
    o_folder = "/mnt/team/msca/pub/archive/evidence-score/gbd2020-process"

    fm = FileManager(i_folder, o_folder)
    pipeline = ContinuousPipeline(fm, PostContinuousProcess)

    for pair in pipeline.pairs:
        print(pair)
        if "metab_bmi_adult" in pair:
            continue
        process = pipeline.get_process(pair)
        process.run()


if __name__ == "__main__":
    main()
