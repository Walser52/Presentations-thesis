<script>
    function showContent(option) {
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
        document.getElementById(option + 'Content').style.display = 'block';

        // Activate the selected nav option
        var selectedNavOption = document.querySelector('.nav-option[data-option="' + option + '"]');
        selectedNavOption.classList.add('active');
    }
</script>
