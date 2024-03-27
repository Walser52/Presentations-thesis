<script src="scripts/jquery-3.6.4.min.js"></script>
<script src="scripts/LeaderLine/leader-line.min.js"></script>




<script>
function toggleBox(header) {
    $(header).next('.dd-box, .ddm-box').slideToggle("slow");
}
</script> 

<script>
    function showContent(clickedElement) {
        // Hide all content divs
        var allContents = document.getElementsByClassName('content');
        for (var i = 0; i < allContents.length; i++) {
            allContents[i].style.display = 'none';
        }

        // Deactivate all nav options
        var allNavOptions = document.getElementsByClassName('nav-option');
        for (var i = 0; i < allNavOptions.length; i++) {
            allNavOptions[i].classList.remove('active');
        }

        // Show the selected content
        var option = clickedElement.getAttribute('data-option');
        document.getElementById(option + 'Content').style.display = 'block';

        // Activate the selected nav option
        clickedElement.classList.add('active');
    }
</script>
